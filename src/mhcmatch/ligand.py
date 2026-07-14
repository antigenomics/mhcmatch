"""Full ligand spans: extend a binding core to the peptide that is actually presented.

Given a 9-mer binding core located in its source protein, :func:`presented_span` returns the most
likely **observed eluted-ligand span** around it -- the peptide a wet lab would synthesise, rather
than the bare core. Three tiers of evidence, weakest last:

* ``observed`` -- a reference ligand in the panel that contains the core *and* occurs in the
  protein. A real eluted span: the gold standard when it exists.
* ``modeled`` -- the highest-scoring feasible span under :class:`SpanModel`, a flank/context model
  fit to mass-spectrometry ligandome data.
* ``fixed`` -- caller-specified flank sizes, clipped at the protein termini (:func:`fixed_span`).

**This is not a cleavage predictor, and not an immunogenicity predictor.** MHC-II peptides are
generated bind-first-trim-later: the groove protects the core while exopeptidases erode the flanks,
so there is no strong sequence-specific endoprotease step to simulate (Paul et al. 2018,
PMID 30127785 -- a dedicated MHC-II cleavage motif reaches AUC 0.767 on ligands and has *zero*
predictive power on CD4 epitopes). What this models is ``P(observed ligand span | source protein)``,
a convolution of protease specificity, HLA-DM editing, binding, stability *and* mass-spectrometry
detection bias. Context/flank models are known to improve *ligand* prediction while **degrading**
CD4 T-cell epitope benchmarks (Reynisson et al. 2020, PMID 32406916). Use this to enumerate and
choose ligands to synthesise or model structurally -- never to rank epitopes by immunogenicity.

For MHC-I the peptide *is* the ligand: there is nothing to extend, so there is no span function.
:func:`processing_score` instead scores an 8-11mer's source-protein context, the shape MHCflurry-2.0
uses for antigen processing (PMID 32711842). Class I and class II are deliberately different
entry points -- a 9-mer class-II core is always <=11 residues and would silently misroute through
any length-based class inference.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from importlib.resources import files

from .diffusion import PROTEOME_AA_FREQ, load_markov1

#: Out-of-protein context position (the span abuts a protein terminus). Modelled, not dropped:
#: a ligand ending exactly at the protein's C-terminus is evidence about where spans end.
PAD = "-"

#: The ligand's own terminal residues -- inside the detected peptide, so subject to MS bias.
LIGAND_KEYS = ("ligN+1", "ligN+2", "ligN+3", "ligC-3", "ligC-2", "ligC-1")
#: The residues flanking the ligand in the source protein -- never in the detected peptide.
FLANK_KEYS = ("flankN-3", "flankN-2", "flankN-1", "flankC+1", "flankC+2", "flankC+3")
#: All 12 context positions: 3 upstream + 3 ligand-N + 3 ligand-C + 3 downstream. This is the
#: NetMHCIIpan ``-context`` window (PMID 30446001); half the signal sits *inside* the ligand.
CTX_KEYS = ("flankN-3", "flankN-2", "flankN-1", "ligN+1", "ligN+2", "ligN+3",
            "ligC-3", "ligC-2", "ligC-1", "flankC+1", "flankC+2", "flankC+3")

CORE_LEN = 9          # MHC-II binding core
_EPS = 1e-6


def _offsets(start, end):
    """Protein index of each context position for the span ``[start, end)``."""
    return {"flankN-3": start - 3, "flankN-2": start - 2, "flankN-1": start - 1,
            "ligN+1": start, "ligN+2": start + 1, "ligN+3": start + 2,
            "ligC-3": end - 3, "ligC-2": end - 2, "ligC-1": end - 1,
            "flankC+1": end, "flankC+2": end + 1, "flankC+3": end + 2}


@dataclass
class Span:
    """A ligand span located in its source protein."""

    peptide: str             # the full ligand, protein[start:end]
    start: int               # 0-based start in the source protein
    end: int                 # 0-based exclusive end
    core: str                # the binding core it was built around
    core_start: int          # 0-based start of the core in the source protein
    source: str              # "observed" | "modeled" | "fixed" -- the evidence tier. ALWAYS check
                             # it: an "observed" span is a lookup, not a prediction.
    score: float = 0.0       # flank/context log-odds (filled for every tier, so tiers compare)
    n_alternatives: int = 0  # other spans within delta of this one. Nested ligand sets are real --
                             # a core usually has several valid spans -- so a lone answer overclaims.
    clipped: tuple = (0, 0)  # (left, right) residues requested but lying outside the protein
    support: int = 0         # reference ligands backing this span ("observed" only)

    @property
    def flanks(self):
        """``(n_left, n_right)`` residues flanking the core within this span."""
        return self.core_start - self.start, self.end - (self.core_start + len(self.core))


@dataclass
class SpanModel:
    """Ligandome-fit flank/context model: ``P(observed ligand span | source protein)``.

    ``ctx`` is allele-agnostic by construction: exopeptidase trimming is a property of the
    proteolytic machinery, not of the groove. That is measured, not assumed -- per-allele context
    PWMs sit within JSD 0.003-0.010 of the pooled one for MHC-II -- and pooling also unlocks the
    ~70% of class-II eluted-ligand records whose restriction is only a placeholder.

    ``lens`` is a **ligand-length** prior, not a core-relative flank-length prior: defining an N-/C-
    flank length requires a binding core, and the allele-agnostic register is tied across >=2 frames
    on ~66% of real ligands, so such a prior would encode a tie-breaking rule rather than biology.
    N/C asymmetry is instead carried by the context positions, which are fit independently per side.

    The span score is a plain log-likelihood, ``log P(L) + context log-odds``, with **no free
    parameters**. A tuned weight on the length prior was tried -- it looked better on the training
    fold and did not transfer (held-out set-recall 0.155 vs 0.158 unweighted, within noise), so it
    was dropped rather than shipped.
    """

    ctx: dict                     # {context key: {residue: frequency}} over CTX_KEYS, incl. PAD
    lens: dict                    # {ligand length: probability}
    padbg: float = 0.02           # marginal rate at which spans abut a terminus -- the PAD null
    background: str = "markov"    # "markov" (order-1 proteome null, conditioned on the preceding
                                  # PROTEIN residue) or "proteome" (order-0). Never the pooled-ligand
                                  # marginal: the training flanks ARE it, so it would be circular.
    _m1: dict = field(default=None, repr=False)

    def __post_init__(self):
        if self.background not in ("markov", "proteome"):
            raise ValueError(
                f"background={self.background!r}: flank log-odds must be scored against a proteome "
                "null. The pooled-ligand marginal is computed FROM these flanks, so it is circular.")
        if self.background == "markov" and self._m1 is None:
            self._m1 = load_markov1()

    def _bg(self, residue, prev):
        """Null probability of ``residue`` given the preceding *protein* residue."""
        if self.background == "markov" and prev and prev in self._m1:
            return self._m1[prev].get(residue) or PROTEOME_AA_FREQ.get(residue, 1e-4)
        return PROTEOME_AA_FREQ.get(residue, 1e-4)

    def context_score(self, protein, start, end, flank_only=False):
        """Log-odds of the context around span ``[start, end)`` vs the proteome null.

        Args:
            protein: the source protein sequence.
            start: 0-based span start.
            end: 0-based exclusive span end.
            flank_only: score only the 6 :data:`FLANK_KEYS`. The 6 ligand-internal positions carry
                the peptide's own anchor signal, so including them partly measures *binding*; the
                flank-only score is the honest processing signal.

        Returns:
            Summed log-odds. Positions outside the protein score against :attr:`padbg`.
        """
        off = _offsets(start, end)
        total = 0.0
        for key in (FLANK_KEYS if flank_only else CTX_KEYS):
            i = off[key]
            f = self.ctx[key]
            if 0 <= i < len(protein):
                r = protein[i]
                prev = protein[i - 1] if i > 0 else ""
                total += math.log((f.get(r, _EPS)) / self._bg(r, prev))
            else:
                total += math.log((f.get(PAD, _EPS)) / self.padbg)
        return total

    def best_span(self, protein, core_start, core_len=CORE_LEN, delta=1.0):
        """Highest-scoring feasible span containing the core, as ``(start, end, score, n_alt)``.

        Only spans that are real substrings of ``protein`` are enumerated, so the result never runs
        off a terminus. ``n_alt`` counts other spans within ``delta`` log-odds of the best -- nested
        sets mean several spans are often legitimately correct.

        The binding term is *identical* for every span sharing this core, so it cancels in the
        argmax and is omitted: ranking is driven purely by the length prior and the flank context.
        (Do not substitute :meth:`AnchorModel.score` here -- it is a max over register frames and so
        grows with peptide length, which would just select the longest span.)
        """
        lens = self.lens
        best, alts = None, []
        for L, pL in lens.items():
            n_flank = L - core_len
            if n_flank < 0:
                continue
            for nl in range(n_flank + 1):
                s = core_start - nl
                e = s + L
                if s < 0 or e > len(protein):
                    continue
                sc = math.log(pL + _EPS) + self.context_score(protein, s, e)
                alts.append(sc)
                if best is None or sc > best[2]:
                    best = (s, e, sc)
        if best is None:                       # core flush against both termini and no L fits
            return core_start, core_start + core_len, 0.0, 0
        n_alt = sum(1 for x in alts if x >= best[2] - delta) - 1
        return best[0], best[1], best[2], n_alt


def _table_path():
    return files("mhcmatch.data").joinpath("ligand_context.tsv")


@lru_cache(maxsize=4)
def load_span_model(cls="mhc2", background="markov"):
    """Load the vendored :class:`SpanModel` for ``cls`` (``"mhc1"`` | ``"mhc2"``).

    Fit from IEDB mass-spectrometry eluted ligands against UniProt reference proteomes; see
    ``src/mhcmatch/data/PROVENANCE.md`` and ``bench/train_spans.py``.
    """
    ctx = {k: {} for k in CTX_KEYS}
    lens, pad_hits, pad_tot = {}, 0.0, 0.0
    with _table_path().open() as fh:
        next(fh)
        for line in fh:
            c, key, b, v = line.rstrip("\n").split("\t")
            if c != cls:
                continue
            if key == "len":
                lens[int(b)] = float(v)
            elif key.startswith("ctx:"):
                k = key[4:]
                ctx[k][b] = float(v)
                if k in FLANK_KEYS and b == PAD:
                    pad_hits += float(v)
                    pad_tot += 1
    if not lens:
        raise ValueError(f"no {cls} rows in {_table_path()}")
    padbg = (pad_hits / pad_tot) if pad_tot else 0.02
    return SpanModel(ctx=ctx, lens=lens, padbg=max(padbg, _EPS), background=background)


def _locate(core, protein, core_start):
    if core_start is None:
        core_start = protein.find(core)
    if core_start < 0 or protein[core_start:core_start + len(core)] != core:
        raise ValueError(f"core {core!r} is not at {core_start} in the given protein")
    return core_start


def observed_spans(core, protein, corpus, core_start=None):
    """Reference ligands that contain ``core`` **and** occur in ``protein``, best first.

    ``corpus`` is any iterable of peptide strings (e.g. ``store._panel['mhc2'].epitopes``). Because
    the caller supplies the protein, ``ligand in protein`` *is* the provenance check -- no source
    accession is needed. Occurrences that do not bracket the core are rejected: a ligand may also
    appear elsewhere in the protein.

    This is a **lookup, not a prediction**. Never fold its hit rate into a prediction metric, and
    never report a training ligand as a novel result -- check :attr:`Span.source`.
    """
    core_start = _locate(core, protein, core_start)
    ce = core_start + len(core)
    hits = {}
    for lig in corpus:
        if core not in lig or lig not in protein:
            continue
        s = -1
        while True:
            s = protein.find(lig, s + 1)
            if s < 0:
                break
            if s <= core_start and ce <= s + len(lig):        # must bracket the core
                hits[(s, s + len(lig))] = hits.get((s, s + len(lig)), 0) + 1
                break
    return [Span(peptide=protein[s:e], start=s, end=e, core=core, core_start=core_start,
                 source="observed", support=n)
            for (s, e), n in sorted(hits.items(), key=lambda kv: (-kv[1], kv[0][1] - kv[0][0]))]


def fixed_span(core, protein, left, right, strict=False, core_start=None):
    """Extend ``core`` by ``left``/``right`` residues, clipped at the protein termini.

    A requested flank that runs off the protein is **reported, not silently shortened**: the
    shortfall lands in :attr:`Span.clipped`.

    Args:
        core: the binding core (must occur in ``protein``).
        protein: the source protein sequence.
        left: residues requested upstream of the core.
        right: residues requested downstream of the core.
        strict: raise instead of clipping when the flank does not fit.
        core_start: 0-based start of ``core``; defaults to its first occurrence.

    Raises:
        ValueError: if ``core`` is not in ``protein``, or ``strict`` and the flank does not fit.
    """
    core_start = _locate(core, protein, core_start)
    ce = core_start + len(core)
    s, e = max(0, core_start - left), min(len(protein), ce + right)
    clipped = (left - (core_start - s), right - (e - ce))
    if strict and any(clipped):
        raise ValueError(f"flank ({left}, {right}) does not fit: clipped by {clipped}")
    return Span(peptide=protein[s:e], start=s, end=e, core=core, core_start=core_start,
                source="fixed", clipped=clipped)


def presented_span(core, protein, model=None, corpus=None, mode="auto", flanks=(3, 3),
                   core_start=None):
    """The most likely presented ligand span around an MHC-II binding ``core``.

    Args:
        core: the 9-mer binding core, located in ``protein``.
        protein: the source protein sequence.
        model: a :class:`SpanModel`; defaults to the vendored MHC-II model.
        corpus: reference ligands for the ``observed`` tier (e.g. panel epitopes). Optional.
        mode: ``"auto"`` (observed -> modeled), or force one of ``"observed"`` | ``"modeled"`` |
            ``"fixed"``. Benchmarks must pass ``"modeled"``: leaving ``observed`` on turns the
            metric into a coverage statistic.
        flanks: ``(left, right)`` for ``mode="fixed"``.
        core_start: 0-based start of ``core``; defaults to its first occurrence. Pass it explicitly
            when the core repeats in the protein.

    Returns:
        A :class:`Span`, or ``None`` when ``mode="observed"`` and no reference ligand contains the
        core -- itself informative: the core has never been eluted.

    Raises:
        ValueError: if ``core`` is not a substring of ``protein``, or is not 9 residues.
    """
    if len(core) != CORE_LEN:
        raise ValueError(f"MHC-II core must be {CORE_LEN} residues, got {len(core)}: {core!r}. "
                         "For an MHC-I peptide use processing_score() -- the peptide is the ligand.")
    core_start = _locate(core, protein, core_start)
    model = model or load_span_model("mhc2")

    if mode in ("auto", "observed") and corpus:
        hits = observed_spans(core, protein, corpus, core_start)
        if hits:
            top = hits[0]
            top.score = model.context_score(protein, top.start, top.end)
            top.n_alternatives = len(hits) - 1
            return top
    if mode == "observed":
        return None
    if mode == "fixed":
        sp = fixed_span(core, protein, flanks[0], flanks[1], core_start=core_start)
        sp.score = model.context_score(protein, sp.start, sp.end)
        return sp

    s, e, sc, n_alt = model.best_span(protein, core_start, len(core))
    return Span(peptide=protein[s:e], start=s, end=e, core=core, core_start=core_start,
                source="modeled", score=sc, n_alternatives=n_alt)


def processing_score(peptide, protein, model=None, flank_only=False, start=None):
    """Source-protein context log-odds of an MHC-I ``peptide`` -- a score, never a span.

    For MHC-I the peptide *is* the ligand, so there is nothing to extend. This scores how
    ligand-like its context looks (the antigen-processing signal MHCflurry-2.0 models,
    PMID 32711842) and composes into ranking, not into an emitted peptide.

    Args:
        peptide: the 8-11mer, located in ``protein``.
        protein: the source protein sequence.
        model: a :class:`SpanModel`; defaults to the vendored MHC-I model.
        flank_only: score only the 6 flanking positions. The ligand's own termini carry its anchor
            signal, so the full 12-position score partly measures binding rather than processing.
        start: 0-based start of ``peptide``; defaults to its first occurrence.

    Raises:
        ValueError: if ``peptide`` is not in ``protein``, or is not 8-11 residues.
    """
    if not 8 <= len(peptide) <= 11:
        raise ValueError(f"MHC-I peptide must be 8-11 residues, got {len(peptide)}: {peptide!r}. "
                         "For an MHC-II core use presented_span().")
    start = _locate(peptide, protein, start)
    model = model or load_span_model("mhc1")
    return model.context_score(protein, start, start + len(peptide), flank_only=flank_only)
