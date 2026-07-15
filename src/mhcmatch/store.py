"""MHC restriction & presentation from a reference epitope panel.

Productionizes the validated reverse-problem method (seqtree ``bench/bench_mhc_guess.py``):
index reference peptides by their anchored *presentation* signature
(:func:`seqtree.layout.presentation_features`), widen the search scope around a query until it
has enough neighbours, then rank presenting alleles by neighbour **vote fraction** and score
**confidence** by a binomial-tail enrichment over the panel background. The vote fraction is the
ranking statistic (robust to panel skew); the enrichment is the non-binder filter.

Significance theory: ``appendix/mhcmatch.tex`` §2-3 (forward per-allele E-value + reverse problem).
"""
from __future__ import annotations

import csv
import gzip
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass

from seqtree import KmerIndex, SearchParams, layout

_CLASS = {"MHCI": "mhc1", "I": "mhc1", "mhc1": "mhc1",
          "MHCII": "mhc2", "II": "mhc2", "mhc2": "mhc2"}
_SPECIES = {"human": "HomoSapiens", "mouse": "MusMusculus"}
_AA = set("ACDEFGHIKLMNPQRSTVWY")
_SCOPES = (0, 1, 2, 3)
_DEFAULT_LENGTHS = {"mhc1": (8, 9, 10, 11), "mhc2": (13, 14, 15, 16, 17, 18)}

PMHC_REPO = "isalgo/pmhc_data"          # public HF dataset holding the reference presentation tables


def fetch_pmhc(tier: str = "full") -> str:
    """Download the pmhc presentation table for ``tier`` from the public HF dataset :data:`PMHC_REPO`
    and return the local cached path.

    Fetches only ``pmhc/pmhc_<tier>.tsv.gz`` (~4-12 MB) — never the other dataset directories — and
    relies on the ``huggingface_hub`` cache, so it downloads once and is instant thereafter. This lets
    a fresh install or a container bootstrap the reference panel with no pre-staged data, which the
    nextflow/Docker deploy depends on. Override with a local ``path=``/``$MHCMATCH_PMHC`` when present.
    """
    from huggingface_hub import hf_hub_download
    return hf_hub_download(repo_id=PMHC_REPO, repo_type="dataset",
                           filename=f"pmhc/pmhc_{tier}.tsv.gz")


_PROTEOME_ALIAS = {"human": "human.fasta.gz", "mouse": "mouse.fasta.gz"}


def fetch_proteome(name: str = "human") -> str:
    """Download a reference proteome FASTA from the public HF dataset :data:`PMHC_REPO` (``proteome/``)
    and return the local cached path.

    ``name`` is ``"human"`` / ``"mouse"`` — the full UniProt proteomes UP000005640 / UP000000589 (for
    source-protein lookup and peptide-flank extraction) — or a pathogen-proteome stem/filename bundled
    in the same dataset (e.g. ``"ecoli_K12_UP000000625"``, for molecular-mimicry sets). Cached by
    ``huggingface_hub``, so it downloads once. Feeds :meth:`mhcmatch.Proteome.from_hf`.
    """
    from huggingface_hub import hf_hub_download
    fname = _PROTEOME_ALIAS.get(name, name if name.endswith(".fasta.gz") else f"{name}.fasta.gz")
    return hf_hub_download(repo_id=PMHC_REPO, repo_type="dataset", filename=f"proteome/{fname}")


def infer_class(peptide: str) -> str:
    """Heuristic class from length: MHC-I if <=11, else MHC-II. Pass ``cls`` to override."""
    return "mhc1" if len(peptide) <= 11 else "mhc2"


@dataclass
class Restriction:
    allele: str
    vote: float        # neighbour vote fraction P(allele | neighbours) -- ranking score
    enrichment: float  # -log10 binomial-tail p vs panel background -- confidence
    n_votes: int
    binder: bool
    anchor_score: float | None = None  # diffused anchor log-odds (set only when diffuse=True)
    rank: float | None = None          # per-allele %rank vs random background (set when calibrated=True)
    p_present: float | None = None     # calibrated presentation probability (set when calibrated=True)
    band: str | None = None            # strong / weak / non-binder from %rank (set when calibrated=True)

    def __iter__(self):
        return iter((self.allele, self.vote, self.enrichment, self.binder))


@dataclass
class Decomposition:
    peptide: str
    tcr_facing: str    # anchors masked with X  (recognition readout)
    presentation: str  # TCR-facing masked with X (anchor readout)
    anchors: tuple     # 0-based anchor indices


def _binom_sf(k, n, p):
    """P(Binomial(n, p) >= k) -- upper tail."""
    if k <= 0:
        return 1.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    return min(1.0, sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1)))


def _bh_cutoff(pvals, alpha):
    """Benjamini-Hochberg p-value cutoff controlling FDR <= ``alpha``: the largest ``p_(r)`` with
    ``p_(r) <= r/m * alpha`` (``0`` if none qualify, so nothing is called)."""
    m = len(pvals)
    if m == 0:
        return 0.0
    cutoff = 0.0
    for r, p in enumerate(sorted(pvals), 1):
        if p <= r / m * alpha:
            cutoff = p
    return cutoff


def _mhc2_register(peptide: str):
    """0-based start of the register-anchored 9-mer core, or None if ``peptide`` is shorter than 9.

    This is the **heuristic** register: a one-pass, allele-agnostic argmax of
    ``seqtree.layout._core_anchor_score`` (leftmost wins ties). It is the register used for
    signatures, ``decompose`` and logos, where no allele is available. The per-allele register the
    model actually scores with is :meth:`mhcmatch.diffusion.AnchorModel.best_register`; on real
    ligands the two often disagree, and both are kept on purpose.
    """
    if len(peptide) < 9:
        return None
    return max(range(len(peptide) - 8),
               key=lambda s: layout._core_anchor_score(peptide[s:s + 9]))


def _mhc2_core_anchors(peptide: str) -> tuple:
    """0-based P1/P4/P6/P9 indices of the register-anchored 9-mer core (one-pass register trick)."""
    s = _mhc2_register(peptide)
    return () if s is None else tuple(s + j for j in (0, 3, 5, 8))


def anchor_indices(peptide: str, cls: str) -> tuple:
    """0-based anchor positions for a peptide: class-I P2/PΩ, class-II core P1/P4/P6/P9."""
    if cls == "mhc2":
        return _mhc2_core_anchors(peptide)
    return tuple(sorted(layout.spec_for(cls).resolve(len(peptide))))


def resolve_anchor_index(peptide: str, cls: str, anchor: int):
    """0-based index of a scoring ``anchor`` in ``peptide`` (or None if out of range).

    MHC-I: ``anchor`` is a 1-based peptide position (negatives count from the C-terminus).
    MHC-II: ``anchor`` is a 1-based position *within the register-anchored 9-mer core* (P1..P9).
    """
    if cls == "mhc2":
        s = _mhc2_register(peptide)
        if s is None:
            return None
        idx = s + (anchor - 1)
        return idx if s <= idx < s + 9 else None
    idx = (anchor - 1) if anchor > 0 else (len(peptide) + anchor)
    return idx if 0 <= idx < len(peptide) else None


class _Panel:
    """One MHC class: presentation-signature KmerIndex + allele bookkeeping."""

    def __init__(self, cls):
        self.cls = cls
        self.epitopes = []
        self.alleles = []
        self.weights = []
        # Unbuilt-panel defaults: a Store used only for decompose() (never loaded via
        # from_records/from_pmhc) must still answer restriction()/alleles() gracefully
        # (empty result) instead of AttributeError-ing on the not-yet-set build() outputs.
        self.index = None
        self.panel = []
        self.freq = {}
        self.allele_to_id = {}

    def add(self, epitope, allele, weight=1.0):
        self.epitopes.append(epitope)
        self.alleles.append(allele)
        self.weights.append(weight)

    def build(self):
        feats = [layout.presentation_features(e, self.cls, register="anchored")
                 for e in self.epitopes]
        self.allele_to_id = {}
        ids = []
        for a in self.alleles:
            self.allele_to_id.setdefault(a, len(self.allele_to_id))
            ids.append(self.allele_to_id[a])
        self.index = KmerIndex.build(feats, alphabet="aa", allele_ids=ids) if feats else None
        counts = Counter(self.alleles)
        total = len(self.alleles) or 1
        self.panel = sorted(counts)
        self.freq = {a: counts[a] / total for a in self.panel}

    def tally(self, query, lo=10, hi=100):
        """Counter(allele -> votes) from the query's anchored-signature neighbours, scope-widened."""
        if self.index is None:
            return None
        feats = layout.presentation_features(query, self.cls, register="anchored")
        cands = []
        for sc in _SCOPES:
            p = SearchParams(max_subs=sc, engine="seqtm")
            cands = [c for c in self.index.seed_and_gather([feats], p, 1, -1, 1)[0]
                     if self.epitopes[c.peptide_id] != query]
            if len(cands) >= lo:
                break
        if not cands:
            return None
        return Counter(self.alleles[c.peptide_id] for c in cands[:hi])


class Store:
    """Searchable reference panel of presented peptides, partitioned by MHC class."""

    def __init__(self):
        self._panel = {"mhc1": _Panel("mhc1"), "mhc2": _Panel("mhc2")}
        self._am = {}  # cls -> AnchorModel (lazy, for diffuse=True)
        self._rc = {}  # cls -> RankCalibrator (lazy, for calibrated=True)

    # -- construction ---------------------------------------------------------
    @classmethod
    def from_records(cls, records):
        """records: dicts with ``epitope``, ``mhc_a`` (or ``mhc``), ``mhc_class``; optional
        ``weight`` (default 1.0) confidence-weights the peptide in anchor-preference estimation."""
        from .pseudoseq import class2_key
        store = cls()
        for r in records:
            c = _CLASS.get(str(r.get("mhc_class", "")).strip())
            ep = str(r.get("epitope", "")).strip().upper()
            allele = str(r.get("mhc_a") or r.get("mhc") or "").strip()
            if c is None or not ep or not allele or not all(x in _AA for x in ep):
                continue
            if c == "mhc2":  # key class II by the alpha-beta pair (locus-aware)
                allele = class2_key(allele, str(r.get("mhc_b") or "").strip())
            store._panel[c].add(ep, allele, float(r.get("weight", 1.0) or 1.0))
        for p in store._panel.values():
            p.build()
        return store

    @classmethod
    def from_pmhc(cls, path=None, tier="full", species=None, classes=("mhc1", "mhc2")):
        """Load the isalgo/pmhc_data TSV(.gz). ``species`` filters the *MHC* species
        (``"human"`` / ``"mouse"``). If ``path`` is None it uses ``$MHCMATCH_PMHC/pmhc_<tier>.tsv.gz``
        when that env var is set, otherwise **bootstraps the table from the public HF dataset** via
        :func:`fetch_pmhc` (downloads only ``pmhc/pmhc_<tier>.tsv.gz``, cached) — so a fresh install or
        a container needs no pre-staged data."""
        if path is None:
            base = os.environ.get("MHCMATCH_PMHC")
            path = os.path.join(base, f"pmhc_{tier}.tsv.gz") if base else fetch_pmhc(tier)
        sp = _SPECIES.get(species) if species else None
        keep = {_CLASS[c] for c in classes}
        csv.field_size_limit(10 ** 7)
        op = gzip.open if str(path).endswith(".gz") else open
        recs = []
        with op(path, "rt") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                c = _CLASS.get(str(row.get("mhc_class", "")).strip())
                if c is None or c not in keep:
                    continue
                if sp and row.get("mhc_species") != sp:
                    continue
                recs.append(row)
        return cls.from_records(recs)

    def __len__(self):
        return sum(len(p.epitopes) for p in self._panel.values())

    def alleles(self, cls):
        return list(self._panel[cls].panel)

    # -- forward problem: restriction / presentation --------------------------
    def _allele_set(self, panel, alleles):
        if alleles == "all":
            return panel.panel
        if isinstance(alleles, str):
            alleles = [alleles]
        return [a for a in alleles if a in panel.freq]

    def _anchor_model(self, cls):
        if cls not in self._am:
            self._am[cls] = self.anchor_model(cls)
        return self._am[cls]

    def _rank_calibrator(self, cls, n=10000, seed=0):
        """Lazy per-allele %rank / P(present) calibrator over a random-peptide background."""
        if cls not in self._rc:
            from .calibrate import RankCalibrator
            panel = self._panel[cls]
            pos = defaultdict(list)
            for ep, a in zip(panel.epitopes, panel.alleles):
                pos[a].append(ep)
            self._rc[cls] = RankCalibrator(self._anchor_model(cls), list(pos), panel.epitopes,
                                           n=n, seed=seed, positives=pos)
        return self._rc[cls]

    def restriction(self, peptide, cls=None, alleles="all", top=10, alpha=0.05, diffuse=False,
                    calibrated=False):
        """Rank presenting alleles for ``peptide`` (vote fraction), flag binders (enrichment).

        ``alleles``: ``"all"``, a single allele, or a list. ``alpha``: per-allele significance for
        the non-binder flag (binder iff binomial-tail p <= alpha and the allele got votes).

        ``calibrated=True`` (implies ``diffuse``) additionally fills each result's ``rank`` (per-allele
        %rank vs a random-peptide background, lower = stronger -- NetMHCpan ``%Rank_EL`` analogue),
        ``p_present``, and qualitative ``band`` (strong/weak/non-binder). The %rank is the
        cross-allele-comparable score; it also re-ranks the results (ascending %rank).

        With ``diffuse=True`` the diffusion-shrunk anchor log-odds
        (:class:`mhcmatch.diffusion.AnchorModel`) **ranks** and the neighbour vote/enrichment
        **gates**: an allele is a binder if it is vote-significant *or* its anchors are plausible
        (``anchor_score > 0``). On held-out (novel) peptides the anchor log-odds is the far better
        ranker---the vote method relies on same-allele signature neighbours, which are sparse for a
        genuinely new peptide, so vote-first ranking buries the true allele; the diffused anchor
        score scores every allele directly and rescues rare ones. Vote breaks ties. Without
        diffusion, vote fraction ranks and the call returns ``[]`` when there are no neighbours.
        """
        peptide = peptide.strip().upper()
        cls = cls or infer_class(peptide)
        panel = self._panel[cls]
        tally = panel.tally(peptide)
        if tally is None and not diffuse:
            return []
        n = sum(tally.values()) if tally else 0
        thr = -math.log10(alpha)
        diffuse = diffuse or calibrated
        am = self._anchor_model(cls) if diffuse else None
        cal = self._rank_calibrator(cls) if calibrated else None
        out = []
        for a in self._allele_set(panel, alleles):
            k = tally.get(a, 0) if tally else 0
            vote = k / n if n else 0.0
            enr = -math.log10(max(_binom_sf(k, n, panel.freq[a]), 1e-300)) if (k and n) else 0.0
            if diffuse:
                s = am.score(peptide, a)
                binder = (enr >= thr and k > 0) or s > 0.0
                r = Restriction(a, vote, enr, k, binder, round(s, 3))
                if cal is not None:
                    from .calibrate import band
                    r.rank = round(cal.percent_rank(a, s), 3)
                    r.p_present = round(cal.p_present(a, s), 4)
                    r.band = band(r.rank)
                out.append(r)
            else:
                out.append(Restriction(a, vote, enr, k, enr >= thr and k > 0))
        if calibrated:
            out.sort(key=lambda r: (r.rank if r.rank is not None else float("inf")))  # ascending %rank
        else:
            out.sort(key=(lambda r: (r.anchor_score, r.vote)) if diffuse
                     else (lambda r: (r.vote, r.enrichment)), reverse=True)
        return out[:top]

    def is_binder(self, peptide, allele, cls=None, alpha=0.05):
        res = self.restriction(peptide, cls=cls, alleles=[allele], top=1, alpha=alpha)
        return bool(res and res[0].binder)

    def is_presented(self, peptide, cls=None, alpha=0.05):
        """Overall presentation: does any panel allele present this peptide?"""
        return any(r.binder for r in self.restriction(peptide, cls=cls, alpha=alpha))

    def scan_protein(self, protein, cls="mhc1", alleles="all", lengths=None, alpha=0.05, top=3,
                     correction=None):
        """Slide all binding-length windows over ``protein`` and return presented peptides.

        Returns ``[(position, peptide, [Restriction, ...]), ...]`` for windows with >=1 binder.

        ``correction`` controls multiple testing over the (window, allele) presentation tests in the
        scan (appendix §5): ``None`` (default) keeps the per-window per-allele ``alpha``;
        ``"bonferroni"`` controls the family-wise error rate (threshold ``alpha/m``); ``"bh"`` controls
        the Benjamini-Hochberg false-discovery rate. ``m`` is the number of voted (window, allele)
        tests. The vote tail p-value is ``10**(-enrichment)``; corrected calls replace the per-window
        binder flag.
        """
        protein = "".join(protein.split()).upper()
        lengths = lengths or _DEFAULT_LENGTHS[cls]
        nA = len(self._allele_set(self._panel[cls], alleles))
        hits = []  # windows that returned candidate alleles, before multiple-testing control
        for L in lengths:
            for i in range(len(protein) - L + 1):
                pep = protein[i:i + L]
                if not all(c in _AA for c in pep):
                    continue
                rs = self.restriction(pep, cls, alleles, top=nA, alpha=alpha)
                if rs:
                    hits.append((i, pep, rs))
        if correction is None:
            return [(i, pep, [r for r in rs if r.binder][:top]) for i, pep, rs in hits
                    if any(r.binder for r in rs)]
        pvals = [10 ** (-r.enrichment) for _, _, rs in hits for r in rs if r.n_votes > 0]
        if correction == "bonferroni":
            cutoff = alpha / len(pvals) if pvals else 0.0
        elif correction == "bh":
            cutoff = _bh_cutoff(pvals, alpha)
        else:
            raise ValueError(f"unknown correction {correction!r} (None|'bonferroni'|'bh')")
        out = []
        for i, pep, rs in hits:
            keep = sorted((r for r in rs if r.n_votes > 0 and 10 ** (-r.enrichment) <= cutoff),
                          key=lambda r: r.enrichment, reverse=True)
            if keep:
                out.append((i, pep, keep[:top]))
        return out

    # -- anchor / TCR-facing split -------------------------------------------
    def decompose(self, peptide, cls=None, allele=None):
        """Split ``peptide`` into anchor and TCR-facing parts, each masked with ``X``.

        ``tcr_facing``: anchors -> X (the recognition readout). ``presentation``: TCR-facing -> X
        (the anchor readout). ``allele`` is accepted for forward-compat (allele-specific learned
        anchors, Phase 1); v0 uses class-default anchor positions.
        """
        peptide = peptide.strip().upper()
        cls = cls or infer_class(peptide)
        anchors = set(anchor_indices(peptide, cls))
        tcr = "".join(layout.MASK if i in anchors else c for i, c in enumerate(peptide))
        present = "".join(c if i in anchors else layout.MASK for i, c in enumerate(peptide))
        return Decomposition(peptide, tcr, present, tuple(sorted(anchors)))

    # -- diffusion-powered forward scorer -------------------------------------
    def anchor_model(self, cls="mhc1", h=2.0, prior_strength=10.0, anchors=None, learn_weights=True,
                     prune_dpi=False, weights="learned", register_em=2, footprint="anchor",
                     rare_max=30, background="ligand"):
        """Anchor-factored presentation model with cross-allele kernel-shrinkage diffusion.

        See :class:`mhcmatch.diffusion.AnchorModel`. The diffusion rescues rare alleles by borrowing
        anchor preferences from groove-similar frequent ones, with a bounded prior strength so a
        large neighbour cannot swamp a rare allele's own peptides. ``register_em`` (MHC-II) runs
        that many best-frame register-EM passes so training and scoring share the same register.
        ``footprint="anchor"`` (default) scores the primary pockets only; ``"core"`` scores the whole
        binding core (MHC-I P1-P5 + PΩ-3..PΩ, MHC-II 9-mer core) -- more discriminative when
        non-anchor positions carry allele-specific signal. ``background="ligand"`` (default) is the
        allele-specificity null; ``"proteome"`` is the presentation null (better for ligand-vs-random
        screening) -- see :data:`mhcmatch.diffusion.PROTEOME_AA_FREQ`.
        """
        from .diffusion import AnchorModel
        return AnchorModel(self, cls=cls, anchors=anchors, h=h, prior_strength=prior_strength,
                           learn_weights=learn_weights, prune_dpi=prune_dpi, weights=weights,
                           register_em=register_em, footprint=footprint, rare_max=rare_max,
                           background=background)

    def affinity_model(self, cls="mhc1"):
        """Quantitative IC50 (nM) + neoantigen amplitude/DAI head (:class:`mhcmatch.PottsAffinity`).

        Loads the vendored Potts weights ``data/affinity_potts_<cls>.npz`` (fields + peptide×pocket
        couplings, fit on measured IEDB IC50). For MHC-II it also builds the register oracle (an
        ``AnchorModel`` with the same ``proteome``/``core`` config used at fit time) so the 9-mer core
        is located consistently. Cached per class. Predict with ``.predict_ic50(peptide, allele)`` and
        the differential ``.amplitude(wt, mut, allele)`` / ``.dai(wt, mut, allele)``.
        """
        from .affinity import PottsAffinity
        cache = self.__dict__.setdefault("_affinity", {})
        if cls not in cache:
            am = self.anchor_model(cls, background="proteome", footprint="core") \
                if cls == "mhc2" else None
            cache[cls] = PottsAffinity(cls, anchor_model=am)
        return cache[cls]

    # -- per-allele anchor preferences (feeds pseudoseq diffusion) ------------
    def anchor_preferences(self, cls, anchor):
        """{allele: Counter(residue)} at a 1-based ``anchor`` position (negative from C-term)."""
        panel = self._panel[cls]
        prefs = defaultdict(Counter)
        for ep, a, w in zip(panel.epitopes, panel.alleles, panel.weights):
            idx = resolve_anchor_index(ep, cls, anchor)
            if idx is not None:
                prefs[a][ep[idx]] += w
        return prefs
