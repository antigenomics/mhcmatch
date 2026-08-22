"""Assembling a polyepitope vaccine cassette: what to refuse, how many units, in what order,
joined by what.

Candidate *selection* — which mutations are worth targeting — is :mod:`mhcmatch.rank` and the
immunogenicity stack behind it. This module is the step after: given ranked candidates with
calibrated probabilities, decide **what to withdraw on safety grounds, how many to carry, how to lay
them out, and what to put between them**. Those are four separate questions with four different
literatures, and none of them is answered by the candidate score.

Selection and assembly are kept apart because the assembly answer depends on the *set*, not the
candidate: whether to carry a 12th epitope depends on what the first eleven already cover, and the
cost of a junction depends on which two units sit either side of it.

**What to refuse — an exclusion, and it runs first.** :func:`screen` withdraws a unit whose own
target gene is transcribed in a tissue that must not be attacked, or one of whose registers coincides
with a self peptide from an unrelated essential-tissue gene. It excludes rather than down-ranks: the
second-best cassette is cheap and myocarditis is not, and capacity spent on a unit that has to be
withdrawn is capacity not spent on a safe one. Two patients died of cardiogenic shock and two of
necrotising leukoencephalopathy in the trials that define this problem; the references, the shape of
each event, and the measurement that chose this screen over the more obvious mimicry-similarity one
are in :func:`screen` and :func:`self_origin_risk`.

**How many — a per-allotype stopping rule, not a constant.**

The clinical numbers (20 for autogene cevumeran across two RNAs, 34 for mRNA-4157, 20 in four pools
for NeoVax) derive from no published objective function, and no trial has delivered N versus 2N
epitopes at matched dose to measure the cost of the extra ones. What *is* established is the shape of
the competition, and it is not a single global budget:

- Competition is **for the antigen-presenting cell, not for MHC**. Suppression required
  co-presentation on the same DC and was reversed by injecting excess pulsed DCs (Kedl et al.,
  *J Exp Med* 2000, PMID 11034600); the mechanism is CD27 cleavage capturing CD70 on the DC
  (Burchill et al., *Eur J Immunol* 2015, PMID 26179759).
- It is **strongest within an allotype and can be net-positive across allotypes**: with large
  pre-existing CTL, a multi-epitope vector "failed to prime efficiently new CTL responses that were
  restricted by the same MHC gene … and vaccine-induced CTL responses restricted by other MHC genes
  were enhanced" (Sherritt et al., *Eur J Immunol* 2000, PMID 10671225).
- The **total response is not a fixed pool**. Deleting four dominant LCMV epitopes gave only "minor
  response increases … and no new epitopes being recognized" (Kotturi et al., *J Immunol* 2008,
  PMID 18641351).

So the expected yield is a sum of **independently saturating per-allotype terms**. With
``p_{a,1} >= p_{a,2} >= ...`` the calibrated probabilities on allotype ``a`` and within-allotype
competition saturating as ``n0 / (n0 + n_a)``::

    E = sum_a  n0 * S_a(n_a) / (n0 + n_a),        S_a(n) = sum_{i<=n} p_{a,i}

    add the next unit on allotype a   <=>   p_{a,n+1} > S_a(n_a) / (n0 + n_a)

:func:`select` implements exactly that line: **keep adding to an allotype while the next candidate's
probability beats that allotype's current expected yield per slot.** Because each allotype saturates
on its own, the marginal value of a crowded allotype's next unit falls below an empty allotype's
first one, so **diversification across allotypes falls out of the arithmetic** instead of being
imposed as a quota. ``n0`` is the one free parameter and it means *per-allotype capacity*; it is not
fitted here, because nothing in the public record fits it. Pass the value you can defend and record
it — :attr:`Selection.n0` carries it into the result so a cassette can always name its own
assumption.

**What between them — scan, do not assume.**

The only linker with causal evidence is ``GPGPG``: unspaced concatenation of four HLA-DR epitopes
created a high-affinity junctional epitope that suppressed the response to *all four*, and inserting
``GPGPG`` restored all four (Livingston et al., *J Immunol* 2002, PMID 12023344). Against that, all
six orderings of three Fel d I regions produced no detectable junctional responses at all (Rogers et
al., *Mol Immunol* 1994, PMID 7521933). Both results are real, which is why this module **measures
each junction against the recipient's own allotypes** rather than picking a linker by reputation.

Two mechanistic constraints worth knowing before choosing a spacer:

- ``AAY`` ends in tyrosine. ERAP1 prefers hydrophobic C-termini and has low affinity for charged
  ones (Chang et al., *PNAS* 2005, PMID 16286653), so a terminal Y genuinely aids processing --
  while also supplying the C-terminal anchor for A\\*01:01, A\\*29:02 and B\\*35:01. It is a
  trade-off, not a mistake, and which way it falls is donor-specific.
- ``KK`` leaves charged residues at the boundary: the mirror image, poor for ERAP1.
- Gly/Pro-rich spacers use residues disfavoured at MHC-I anchor positions and abundant in the
  C-terminal regions from which ligands are cleaved (Martin-Galiano & Lopez, *PLoS One* 2019,
  PMID 30645615), so they sit in the permissive zone.

:data:`SPACERS` therefore leads with ``None``. A cassette whose junctions are already clean should
carry no spacer at all — every inserted residue is sequence that has to be translated and could
itself form a binder.

**Order** is chosen the way ``pVACvector`` chooses it (Hundal et al., *Cancer Immunol Res* 2020,
PMID 31907209): score every register spanning every junction, treat the units as nodes of a complete
graph whose edge cost is the *strongest* predicted binder at that junction, and find a cheap open
path. This module uses a deterministic greedy path plus bounded 2-opt rather than simulated
annealing — no RNG, so a cassette is reproducible from its inputs.

**Scoring is injected, never imported.** :func:`order` and :func:`scan_junctions` take a
``binder`` callable, so the layout logic is testable without a :class:`~mhcmatch.store.Store`, a
panel or any download, and a caller who wants a different presentation model just passes it::

    >>> from mhcmatch import vector
    >>> units = [vector.Unit("AAAAAAAAAKAAAAAAAAAAAAAAAAA", 9, "GENE1", "HLA-A*02:01", 0.30),
    ...          vector.Unit("CCCCCCCCCRCCCCCCCCCCCCCCCCC", 9, "GENE2", "HLA-A*02:01", 0.12)]
    >>> sel = vector.select(units, n0=8)
    >>> [u.gene for u in sel.units]
    ['GENE1', 'GENE2']
    >>> cas = vector.order(sel.units, binder=lambda peps, alleles: [0.0] * len(peps))
    >>> cas.spacer is None and len(cas.units) == 2
    True
    >>> cas.sequence == units[0].peptide + units[1].peptide
    True

A ``binder`` returns one number per peptide, **higher meaning a stronger predicted binder**;
``-log10(%rank)`` is the natural choice and is what :func:`store_binder` builds.

From the command line, where the whole pipeline is one call::

    mhcmatch vector --candidates units.tsv --n0 8 --screen --fasta cassette.fasta

``units.tsv`` carries ``peptide``, ``gene``, ``allele``, ``p`` and optionally ``mutation_index`` —
and ``peptide`` is the **long window**, not ``rank``'s minimal epitope. ``--screen`` is opt-in
because it builds a whole-proteome index; without it no safety check runs at all and the cassette
carries whatever it was handed. ``mhcmatch deslip <cds> --fix out.fasta`` is the
:func:`slippery_sites` half, which takes nucleotides rather than peptides and so is its own command.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

from .predict import NOVEL_PRODUCTS, TRACT_PRODUCTS

#: Spacers tried in order, ``None`` first. A clean junction needs no spacer, and pVACvector's
#: default list is tried only when one is needed. See the module docstring for why ``AAY`` is a
#: trade-off rather than a default.
SPACERS: tuple = (None, "GPGPG", "GGS", "AAA", "AAY", "GPGPGPG", "HHAA", "AAL")

#: MHC-I register lengths scanned across a junction. Class II cores are 9-mers read out of a longer
#: span, so a class-II junction scan wants :data:`MHC2_JUNCTION_LENGTHS` instead.
JUNCTION_LENGTHS: tuple = (8, 9, 10, 11)

#: Class-II junction registers. The bound core is 9 residues but the presented span is longer, so
#: the scan needs the window a core could be read from.
MHC2_JUNCTION_LENGTHS: tuple = (12, 13, 14, 15)

#: The 20 standard amino acids, in the order the substitution neighbourhood is enumerated.
_AA: str = "ACDEFGHIKLMNPQRSTVWY"


def _neighbours(peptide: str) -> list:
    """Every sequence exactly one substitution from ``peptide`` (``19 * len(peptide)`` of them).

    The peptide itself is excluded, so a hit on a neighbour is a genuine ``d = 1`` coincidence and
    never the exact match :meth:`~mhcmatch.Proteome.find_exact_sources` would already have found.
    """
    return [peptide[:i] + a + peptide[i + 1:]
            for i in range(len(peptide)) for a in _AA if a != peptide[i]]


def flank_identity(a: str, ai: int, b: str, bi: int, length: int, k: int = 10) -> float:
    """Fraction of matching residues in the ``k`` positions either side of a shared register.

    ``a``/``b`` are the two contexts the register was found in (a unit's 27-mer and a reference
    protein), ``ai``/``bi`` the register's 0-based start in each, ``length`` its length. Positions
    truncated by either sequence's end are not compared, and a register flush against both ends
    scores 0.0 — nothing was compared, so nothing supports homology.

    **This is what separates a paralog from a coincidence, and the separation is the whole reason
    the report tier is usable.** Two proteins from one locus share their flanks as well as the
    register, so a ``d = 1`` match between them is descent, not mimicry, and the T cell that sees one
    is a T cell tolerance already had to deal with. A match whose flanks are unrelated is an
    independent occurrence: the same nine residues arrived twice by chance, and *that* is the object
    the titin and MAGE-A12 deaths are drawn from. Measured on 178 validated immunogenic neoantigens,
    the cut at 0.5 removes **156 of 230** different-gene ``d = 1`` hits at L=9 — 130 of them at 90%
    or better, which is one locus under two symbols rather than a coincidence at all.

    The comparison is bounded by whatever context the caller gives, and a 27-mer unit gives ±9-10
    residues. That is enough to separate **loci**, not superfamilies: NRAS → KRAS scores 0.23 and is
    reported, which is the wanted behaviour — a T cell raised on an NRAS Q61 neoantigen that
    cross-reacts to wild-type KRAS is a real on-target/off-tumour concern, and KRAS is transcribed
    everywhere.
    """
    lo = min(ai, bi, k)
    hi = min(len(a) - ai - length, len(b) - bi - length, k)
    n = lo + hi
    if not n:
        return 0.0
    same = sum(x == y for x, y in zip(a[ai - lo:ai], b[bi - lo:bi]))
    same += sum(x == y for x, y in zip(a[ai + length:ai + length + hi],
                                       b[bi + length:bi + length + hi]))
    return same / n

#: GTEx ``SMTSD`` tissue **prefixes** whose destruction is not survivable or not repairable, so a
#: candidate whose self-mimic is transcribed there is excluded rather than ranked down. Matched by
#: prefix because GTEx splits organs into regions -- ``Brain`` alone covers twelve of the 53 tissue
#: names, ``Heart`` two.
#:
#: The list is a clinical judgement, not a fitted parameter, and both fatal precedents behind it are
#: in :func:`screen`'s docstring. Override it: a cassette for a patient who has already lost an organ
#: is a different calculation, and so is one for a tissue the protocol accepts damaging.
#: The nine organs whose damage is not survivable or not recoverable. **Both naming schemes.**
#: The expression table carries 123 distinct context names from two sources -- GTEx-style
#: (``"Brain - Caudate (basal ganglia)"``, ``"Muscle - Skeletal"``) and HPA-style lowercase
#: (``"basal ganglia"``, ``"skeletal muscle"``) -- and a nine-entry Title-Case tuple matched by
#: ``str.startswith`` sees only **22 of the 123**. Thirteen essential organs were invisible:
#: heart muscle, kidney, liver, lung, adrenal gland, pituitary gland, cerebellum, cerebral cortex,
#: midbrain, hippocampal formation, spinal cord, skeletal muscle, smooth muscle.
#:
#: Measured cost of that omission, on the shipped table: of the 7,527 genes reaching 50 TPM in an
#: essential tissue, **1,517 (20.2 %) were invisible to the screen** -- among them CEACAM5 (4.65 seen
#: vs 28.50 actual; Parkhurst 2011 colitis, 3 of 3 patients), CDH13 (7.09 vs 64.20) and albumin
#: (26,217 vs 198,524). This is a false negative in the fatal direction, so the list is not a scope
#: choice: it is the same nine organs, spelled the way the data spells them.
ESSENTIAL_TISSUES: tuple = ("Heart", "Brain", "Nerve", "Lung", "Liver", "Kidney",
                            "Adrenal Gland", "Pituitary", "Muscle - Skeletal",
                            # the same organs under the table's lowercase naming
                            "heart muscle", "kidney", "liver", "lung", "adrenal gland",
                            "pituitary gland", "skeletal muscle", "smooth muscle",
                            "cerebellum", "cerebral cortex", "midbrain", "basal ganglia",
                            "hippocampal formation", "hypothalamus", "amygdala",
                            "spinal cord", "choroid plexus", "retina")


@dataclass(frozen=True)
class Unit:
    """One vaccine unit: the long peptide carrying one mutation, plus what it is worth.

    ``peptide`` is placed into the cassette verbatim, so build it with :func:`unit` rather than
    slicing by hand — the mutation has to sit far enough from both ends that every register
    containing it is generated (see :func:`unit`).

    ``allele`` is the restriction the unit is *credited to* for the per-allotype budget in
    :func:`select`. A long peptide usually presents on several allotypes; credit it to the one whose
    probability ``p`` refers to, and carry a second :class:`Unit` for the same mutation only if a
    different allotype genuinely needs a different window.

    ``p`` is a calibrated probability of eliciting a detectable response, on the operating prior the
    caller intends to deploy at — not a corpus-prevalence posterior, which overstates it. Nothing
    here checks that; it is the caller's contract.
    """

    peptide: str
    mutation_index: int
    gene: str
    allele: str
    p: float
    cls: str = "mhc1"
    #: What kind of variant produced the neoepitope. ``"missense"`` is a single substitution against
    #: a self protein; anything else -- ``frameshift``, ``fusion``, ``splice``, ``retained_intron``,
    #: ``ORF``, ``editing`` -- is **non-conventional** and is charged to its own arm by
    #: :func:`mhcmatch.portfolio.compose`. The distinction earns a quota of its own because a
    #: non-conventional product is foreign over a stretch rather than at one position, so it fails
    #: differently from a missense: whatever makes the missense arm miss (a wrong wild type, a
    #: tolerised residue) does not make this arm miss.
    kind: str = "missense"

    def __post_init__(self):
        if not 0 <= self.mutation_index < len(self.peptide):
            raise ValueError(f"mutation_index {self.mutation_index} outside "
                             f"peptide of length {len(self.peptide)}")
        if not 0.0 <= self.p <= 1.0:
            raise ValueError(f"p must be a probability, got {self.p}")


@dataclass
class Selection:
    """What :func:`select` kept, what it dropped, and the rule that decided.

    ``trace`` is one row per *considered* candidate in the order the rule saw it, carrying the
    threshold it was compared against. A cassette that cannot explain why its 12th unit is in and the
    13th is out is not auditable, and the threshold is cheap to record.
    """

    units: list = field(default_factory=list)
    dropped: list = field(default_factory=list)
    n0: float = 0.0
    trace: list = field(default_factory=list)
    keys: list = field(default_factory=list)

    def _grouped(self) -> dict:
        """``{block key: [p, ...]}`` for the blocks the rule actually charged against.

        Falls back to allotype when ``keys`` is absent, so a :class:`Selection` unpickled from
        before the ``block`` parameter existed still reports the same numbers it always did.
        """
        keys = self.keys if len(self.keys) == len(self.units) else [u.allele for u in self.units]
        by: dict = {}
        for key, u in zip(keys, self.units):
            by.setdefault(key, []).append(u.p)
        return by

    @property
    def expected_yield(self) -> float:
        """``sum_b n0 * S_b / (n0 + n_b)`` — expected responses under the saturation model.

        Grouped by whatever :func:`select` blocked on, not unconditionally by allotype: the yield
        has to be computed against the same partition the stopping rule spent its budget on, or it
        describes a cassette that was never built.
        """
        return sum(self.n0 * sum(ps) / (self.n0 + len(ps)) for ps in self._grouped().values())

    def per_allele(self) -> dict:
        """``{allele: (n_units, summed p, saturated yield)}`` — where the budget actually went."""
        by: dict = {}
        for u in self.units:
            by.setdefault(u.allele, []).append(u.p)
        return {a: (len(ps), sum(ps), self.n0 * sum(ps) / (self.n0 + len(ps)))
                for a, ps in sorted(by.items())}

    def per_block(self) -> dict:
        """``{block key: (n_units, summed p, saturated yield)}`` for the rule's own partition."""
        return {k: (len(ps), sum(ps), self.n0 * sum(ps) / (self.n0 + len(ps)))
                for k, ps in sorted(self._grouped().items(), key=lambda kv: str(kv[0]))}


@dataclass
class Cassette:
    """An ordered, spaced cassette and the junction evidence behind its layout.

    ``sequence`` is the epitope cassette **only** — no start codon, no stop, no leader and no
    trafficking domain. Those belong to the vector backbone, and a cassette that silently included
    them could not be cloned into one that already had them.
    """

    units: list
    spacer: str | None
    sequence: str
    boundaries: list = field(default_factory=list)
    junctions: list = field(default_factory=list)
    cost: float = 0.0

    @property
    def worst_junction(self) -> float:
        """Strongest predicted junctional binder anywhere in the cassette."""
        return max((j["score"] for j in self.junctions), default=float("-inf"))


def unit(context: str, mutation_offset: int, length: int = 27,
         **kw) -> Unit:
    """Centre ``mutation_offset`` of ``context`` in a ``length``-mer window, and wrap it as a
    :class:`Unit`.

    ``length`` defaults to **27** with the mutation at position 14 — the configuration BioNTech's
    backbone carries (Kreiter et al., *Nature* 2015, PMID 25901682) and the one that guarantees every
    8-to-14-residue register containing the mutation is present, so all lengths and all allotypes are
    covered by a single unit and duplication buys nothing.

    Centring is also why a *minimal* epitope should never be a unit. A minimal peptide can load
    directly onto MHC-I of any cell, including non-professional APCs with no costimulation, and does:
    injected alone it "transiently activated CD8+ effector T cells, which eventually failed to undergo
    secondary expansion or to kill target cells", while simply extending it to 30 residues restored
    both, "independent of T cell help, because the longer CTL peptide was predominantly presented in
    the locally inflamed draining lymph node" (Bijker et al., *J Immunol* 2007, PMID 17911588). Short
    units are not merely less efficient, they are the tolerising configuration.

    Near a protein terminus the window is clamped and the mutation sits off-centre; that is
    unavoidable and :attr:`Unit.mutation_index` records where it actually landed.
    """
    if not 0 <= mutation_offset < len(context):
        raise ValueError(f"mutation_offset {mutation_offset} outside context "
                         f"of length {len(context)}")
    half = (length - 1) // 2
    start = max(0, min(mutation_offset - half, len(context) - length))
    start = max(start, 0)
    window = context[start:start + length]
    return Unit(peptide=window, mutation_index=mutation_offset - start, **kw)


def units_from_context(rows, records, length: int = 27, cls: str = "mhc1") -> list:
    """``[Unit]`` from ranked **minimal epitopes** plus the window FASTA they were called on.

    This is the join between :func:`mhcmatch.rank.rank_fasta` and this module. ``rank`` emits
    minimal epitopes and a score; a unit is the long window around the mutation, and where that
    mutation sits is in the FASTA header rather than in ``rank``'s output -- so neither side alone
    can build one. ``rows`` are dicts carrying ``peptide`` (the minimal epitope), ``gene``,
    ``allele`` and ``p``; ``records`` is :func:`mhcmatch.predict.parse_fasta`'s output for the very
    FASTA ``rank`` was pointed at.

    **One unit per variant, not per epitope.** Twenty registers of one mutation are twenty rows in
    ``rank`` and one thing to put in a cassette, and :func:`select` spends capacity per unit -- so
    rows are grouped by their source window and the group's best-scoring row supplies the allotype
    and the score. Anything whose epitope matches no window is returned to the caller's attention by
    being absent; counts are the caller's to report.

    **All four header families are admitted, not only** ``Somatic:``. Each carries its novel residue
    or span somewhere different, and :func:`_centred_context` reads all four; skipping the other
    three discarded 317 of the 489 non-missense records in one real cohort and left the
    ``nonconventional`` quota arm with nothing to fill itself from -- which is the one arm whose
    whole point is that it fails differently from the missense arm.
    """
    from .predict import parse_variant_header, variant_product

    contexts = []
    for header, seq in records:
        var = parse_variant_header(header)
        got = _centred_context(var, seq)
        if got is not None:
            contexts.append((header, var, got[0], got[1]))

    # Indexed, not scanned. `pep in context` over every context per row is O(rows x windows)
    # substring tests -- 4,606 windows against tens of thousands of ranked rows. A context that
    # contains `pep` contains its first K residues, so the K-mer bucket is a superset of the
    # matches and is the only thing that has to be tested; the gene index serves the short peptides
    # a K-mer key cannot be built for.
    K = 8
    by_kmer: dict = {}
    by_gene: dict = {}
    for i, (_h, var, ctx, _o) in enumerate(contexts):
        by_gene.setdefault(var.get("gene_name") or "", []).append(i)
        for j in range(len(ctx) - K + 1):
            by_kmer.setdefault(ctx[j:j + K], []).append(i)

    best: dict = {}
    for r in rows:
        pep = str(r.get("peptide", "")).strip().upper()
        if not pep:
            continue
        gene = str(r.get("gene", "")).strip()
        if len(pep) >= K:
            cand = by_kmer.get(pep[:K], ())
        else:
            cand = by_gene.get(gene) or range(len(contexts))
        hit = None
        for i in cand:
            c = contexts[i]
            if pep not in c[2]:
                continue
            if not gene or c[1].get("gene_name") == gene:
                hit = c
                break
            if hit is None:                      # a match on the wrong gene, kept only as fallback
                hit = c
        if hit is None:
            continue
        p = float(r.get("p", 0.0) or 0.0)
        if hit[0] not in best or p > best[hit[0]][0]:
            best[hit[0]] = (p, r, hit)

    out = []
    for p, r, (_h, var, context, offset) in best.values():
        out.append(unit(context, offset, length=length,
                        gene=var["gene_name"] or str(r.get("gene", "")).strip(),
                        allele=str(r.get("allele", "")).strip(), p=p,
                        cls=str(r.get("cls", "") or cls),
                        # NOT `var["type"]`: that is provenance (`Somatic`), and reading it here
                        # charged every candidate to the non-conventional arm. `rank` emits the
                        # product class in `variant_type`; where the caller's rows carry neither,
                        # the header's own product is read rather than defaulting to `missense` --
                        # a fusion silently charged to the missense arm is the same bug in reverse.
                        kind=str(r.get("kind", "") or r.get("variant_type", "")
                                 or variant_product(var) or "missense")))
    return out


def _centred_context(var: dict, seq: str):
    """``(context, offset)`` for one parsed window header, or ``None`` if it carries no centre.

    ``context`` is the sequence a unit is cut from and ``offset`` the position in it that the unit
    must be built around. The four families mark that position four different ways, all verified
    against real pipeline records:

    - ``Somatic:`` -- the mutated residue is parenthesised in ``mut_window``, so the marker's index
      is the offset and the context is the window with the parens stripped.
    - ``Fusion:`` / ``CNV:`` -- the ``windows`` field is pipe-separated. **Three parts**
      (``LEFT|X|RIGHT``) marks an altered residue at ``len(LEFT)``; **two parts** (``LEFT|RIGHT``)
      marks the junction itself, also at ``len(LEFT)``. A fusion read through a stop codon carries a
      ``*``, and everything past it is not translated, so the context is truncated there.
    - ``Isoform:`` -- no residue is altered; the whole novel span is. The header's trailing ``span``
      field ends in a ``d-e`` pair giving that span **in the emitted sequence**, and the unit is
      centred on it.

    The emitted FASTA sequence is a *sub*-window of the ``windows`` field for every ``Fusion:`` and
    ``CNV:`` record measured (79 of 79), so the pipe-separated field is the context and the emitted
    sequence is not -- taking the latter would put the offset in the wrong coordinate frame.
    """
    from .predict import _strip_marker

    kind = str(var.get("type") or "")
    if kind == "Somatic":
        marked = var.get("mut_window") or ""
        return (_strip_marker(marked), marked.index("(")) if "(" in marked else None
    if kind in ("Fusion", "CNV"):
        parts = (var.get("windows") or "").split("|")
        if len(parts) not in (2, 3):
            return None
        ctx, off = "".join(parts).split("*")[0], len(parts[0])
        return (ctx, off) if 0 <= off < len(ctx) else None
    if kind == "Isoform":
        import re
        nums = [n for n in re.split(r"[:-]", var.get("span") or "") if n.isdigit()]
        if len(nums) < 2 or not seq:
            return None
        mid = (int(nums[-2]) + int(nums[-1])) // 2
        return seq, max(0, min(mid, len(seq) - 1))
    return None


def junction_windows(left: str, right: str, spacer: str | None = None,
                     lengths=JUNCTION_LENGTHS) -> list:
    """Every ``lengths``-mer that spans the ``left``/``right`` boundary, as ``(peptide, offset)``.

    A window *spans* the boundary when it contains at least one residue from each side, so windows
    lying wholly inside either unit are excluded — those are the intended epitopes, not artefacts of
    concatenation. With a spacer, a window containing only spacer plus one side still counts: it did
    not exist in either unit and did not exist in any genome.

    ``offset`` is the window's start relative to the start of ``left``, so a caller can point at the
    exact residues to change.
    """
    sp = spacer or ""
    joined = left + sp + right
    lo, hi = len(left), len(left) + len(sp)      # [lo, hi) is the spacer
    out = []
    for L in lengths:
        for i in range(0, len(joined) - L + 1):
            j = i + L
            if i < lo and j > hi:                # touches both units
                out.append((joined[i:j], i))
            elif sp and i < hi and j > lo:       # touches the spacer and at least one unit
                out.append((joined[i:j], i))
    return out


def scan_junctions(units, binder, spacer: str | None = None,
                   lengths=JUNCTION_LENGTHS, alleles=None,
                   binder_threshold: float | None = None) -> list:
    """Score every junction of ``units`` laid out in the given order.

    Returns one dict per junction: ``left``, ``right`` (unit indices), ``score`` (the **strongest**
    predicted binder spanning it), ``peptide`` and ``offset`` of that worst window, and ``n_windows``.

    The strongest rather than the mean, because a junction is a hazard if it forms *one* good binder;
    averaging it against the many bad windows around it hides exactly the case being looked for. This
    is pVACvector's "lowest binding score" convention read on a scale where higher is a better
    binder.
    """
    out = []
    for k in range(len(units) - 1):
        wins = junction_windows(units[k].peptide, units[k + 1].peptide, spacer, lengths)
        if not wins:
            out.append({"left": k, "right": k + 1, "score": float("-inf"),
                        "peptide": None, "offset": None, "n_windows": 0})
            continue
        peps = [w for w, _ in wins]
        scores = list(binder(peps, alleles))
        if len(scores) != len(peps):
            raise ValueError(f"binder returned {len(scores)} scores for {len(peps)} peptides")
        b = max(range(len(peps)), key=lambda i: scores[i])
        n_over = (sum(1 for v in scores if v >= binder_threshold)
                  if binder_threshold is not None else None)
        out.append({"left": k, "right": k + 1, "score": float(scores[b]),
                    "peptide": peps[b], "offset": wins[b][1], "n_windows": len(peps),
                    "n_over": n_over})
    return out


def screen(units, risk, lengths=JUNCTION_LENGTHS, notes=None) -> tuple:
    """``(kept, rejected)`` — drop units carrying a register that mimics an essential-tissue self
    peptide. **Run this before :func:`select`**: capacity spent on a unit that has to be withdrawn is
    capacity not spent on a safe one.

    The precedent is not hypothetical and it is not a binding-prediction failure. An
    affinity-enhanced TCR against the HLA-A\\*01:01-restricted MAGE-A3 epitope ``EVDPIGHLY`` killed
    the first two patients infused, by cardiogenic shock within days; autopsy found T-cell infiltrate
    and myocardial damage with **no MAGE-A3 expressed in heart at all**, and the off-target was
    ``ESDPIVAQY`` from titin (Linette et al., *Blood* 2013;122(6):863-71, PMID 23770775; Cameron et
    al., *Sci Transl Med* 2013;5(197):197ra103, PMID 23926201). Separately, a TCR recognising
    MAGE-A3/A9/A12 caused necrotising leukoencephalopathy and two deaths, because MAGE-A12 turned out
    to be transcribed in human brain (Morgan et al., *J Immunother* 2013;36(2):133-51, PMID 23377668).

    Read the two together and they give this function its shape:

    - **The hazard is in the source gene's tissue, not the candidate's score.** Both events were
      invisible to binding prediction and visible in expression: MAGE-A12 is transcribed in brain, and
      titin in heart. So the check joins a peptide to a *protein* to a *tissue*.
    - **Two different questions, and the unit answers both.** Is the unit's own target gene
      transcribed somewhere it must not be attacked (MAGE-A12)? And does any register of it coincide
      with a self peptide from an **unrelated** essential-tissue gene?
    - **Exclusion, not down-ranking.** The second-best cassette is cheap; myocarditis is not.

    **The unit's own gene has to be excluded from the register test, or the screen rejects
    everything.** A 27-mer is native context by design: its flanking registers *are* self peptides
    from its own parent protein, and the mutated register sits one substitution from that protein's
    wild type. Screened naively at any useful radius, every unit of every cassette fires. Those
    matches are also the ones tolerance already covers — the flanks are presented in normal tissue
    daily. What is not covered is a register that coincides with a *different* protein, which is why
    ``risk`` is handed the unit and not just its registers.

    ``risk(unit, registers) -> [reason, ...]`` returns **zero or more** reasons, empty meaning safe.
    Each reason is a dict; a ``"register"`` key naming which register triggered it is carried into the
    record, and its absence means the reason is unit-level. :func:`self_origin_risk` builds one from
    the human proteome and reference expression. Injected for the same reason ``binder`` is — the
    policy here is testable with no panel, no proteome and no download, and a site with its own
    toxicity list substitutes it wholesale.

    ``rejected`` is ``[(unit, register, reason)]``, ``register`` being ``None`` for a unit-level
    reason. A withdrawn candidate has to say what withdrew it: "the screen dropped 3 of 40" is not a
    safety argument, and the reason is what a clinician overrides or accepts.

    **A reason may decline to withdraw.** ``r["veto"] = False`` marks a finding that is recorded but
    does not exclude -- the graded mode of :func:`self_origin_risk`, where a hit below ``veto_tpm``
    is a *cost* to composition rather than a refusal. The key's absence means ``True``, so a risk
    callable that never sets it behaves exactly as it always did. Pass ``notes=[]`` to collect those
    non-vetoing findings; they arrive in the same ``(unit, register, reason)`` shape as ``rejected``
    and are what :func:`offtarget_cost` reads.

    **One batch query for the whole candidate list, not one per unit.** Where ``risk`` exposes a
    ``prepare(registers)`` (:func:`self_origin_risk` does), it is handed the deduplicated registers
    of *every* unit before any unit is judged. A 27-mer carries ~70 registers and a real cohort
    carries thousands of units, so the per-unit call pattern made ~19,000 proteome queries where 1
    suffices, and units share registers heavily.

    A junction can manufacture a self-mimic just as it can manufacture a binder, but that check needs
    the layout, so it belongs after :func:`order` rather than here.
    """
    per_unit = [sorted({u.peptide[i:i + L] for L in lengths
                        for i in range(len(u.peptide) - L + 1)}) for u in units]
    prepare = getattr(risk, "prepare", None)
    if prepare is not None:
        prepare(sorted({r for regs in per_unit for r in regs}))

    kept, rejected = [], []
    for u, regs in zip(units, per_unit):
        reasons = list(risk(u, regs)) if regs else []
        vetoes = [r for r in reasons if r.get("veto", True)]
        if vetoes:
            rejected.extend((u, r.get("register"), r) for r in vetoes)
        else:
            kept.append(u)
        if notes is not None:
            notes.extend((u, r.get("register"), r) for r in reasons if not r.get("veto", True))
    return kept, rejected


def offtarget_cost(findings) -> dict:
    """``{unit: cost}`` from :func:`screen`'s ``notes`` (or its ``rejected``) -- the size of a unit's
    **off-target fingerprint**, being the number of distinct ``(clause, gene)`` pairs it reaches.

    Distinct genes, not reasons: a gene is transcribed in many tissues and :func:`self_origin_risk`
    reports one finding per tissue, so counting reasons would charge a unit for the breadth of GTEx
    rather than the breadth of its off-targets. Units absent from ``findings`` are absent from the
    dict; read it with ``.get(u, 0.0)``.

    This is the number :func:`mhcmatch.portfolio.compose` subtracts under ``weight_cost``. It is a
    count and not a probability on purpose -- there is no calibration behind "how much worse is two
    off-target genes than one", so the weight is the caller's to set and to record.
    """
    by: dict = {}
    for u, _reg, why in findings:
        by.setdefault(u, set()).add((why.get("clause", ""), why.get("gene", "")))
    return {u: float(len(s)) for u, s in by.items()}


def presented(findings, binder, threshold: float = -1.47712, alleles=None) -> list:
    """Keep the near-identity ``findings`` whose off-target variant is **actually presented**.

    A ``d = 1`` coincidence is only a hazard if a T cell can see it, and seeing it means the
    *off-target's own* sequence — not the unit's register — is presented on the allotype the unit was
    selected for. A variant that no allotype presents is a sequence coincidence and nothing more, so
    it is dropped from the fingerprint rather than reported as a safety consideration. Findings with
    no ``"variant"`` key (clauses 1 and 2, and any sub-veto finding under ``graded``) pass through
    untouched: presentation is not what they are about.

    ``binder(peptides, alleles) -> [score]`` is :func:`store_binder`'s contract, ``score`` being
    ``-log10(%rank)`` so higher is stronger. Every finding on one allotype goes in **one** call — the
    alternative is a :meth:`~mhcmatch.store.Store.restriction` per finding, and a cassette's report
    tier carries thousands. ``alleles`` overrides the per-unit allotype for callers whose units carry
    none.

    **``threshold`` defaults to ``-log10(30)``, and the conventional 2% rank would be wrong here.**
    This gate is a safety read-out, so the expensive error is missing a hazard, and the cut belongs
    where the *positives* are rather than at a number borrowed from a different scorer. Measured on
    the 176 assayed immunogenic neoantigen/allotype pairs in ``isalgo/pmhc_data``, scored by this
    scorer on their own allotype, the median sits at **0.69% rank** and the 5th percentile at
    **14.3%**:

    ==========  ====================================  ==========================
    %rank cut   assayed immunogenic peptides kept     units carrying a report
    ==========  ====================================  ==========================
    none        100.0%                                27 of 174 (15.5%)
    30          **97.2%**                             **14 of 174 (8.0%)**
    20          96.0%                                 12 of 174 (6.9%)
    15          94.9%                                  9 of 174 (5.2%)
    5           88.6%                                  4 of 174 (2.3%)
    2           70.5%                                  0 of 174 (0.0%)
    ==========  ====================================  ==========================

    At 2% the gate discards **three in ten** genuinely immunogenic peptides, which on a safety
    question is the error that costs something. 30% still halves the tier, 27 units to 14.

    Each finding gains a ``"variant_binder"`` key with its score, kept even when it fails, because
    "we looked and it is not presented" is the part of a safety argument that a bare absence cannot
    make.
    """
    by: dict = {}
    for i, (u, _reg, why) in enumerate(findings):
        if "variant" in why:
            by.setdefault(tuple(alleles) if alleles else (getattr(u, "allele", None),), []).append(i)
    ok = set()
    for allele, idx in by.items():
        scores = binder([findings[i][2]["variant"] for i in idx],
                        [a for a in allele if a] or None)
        for i, score in zip(idx, scores):
            findings[i][2]["variant_binder"] = score
            if score >= threshold:
                ok.add(i)
    return [f for i, f in enumerate(findings) if "variant" not in f[2] or i in ok]


def self_origin_risk(proteome, symbols, tissues=ESSENTIAL_TISSUES, min_tpm: float = 0.25,
                     max_subs: int = 0, *, novel_kinds=NOVEL_PRODUCTS,
                     veto_tpm: float = 5.0, graded: bool = False,
                     report_subs: int = 0, report_identity: float = 0.5,
                     report_flank: int = 10, report_min_length: int = 9):
    """A ``risk`` callable for :func:`screen`: **near-exact self origin, joined to tissue.**

    A register is risky when :meth:`mhcmatch.Proteome.find_source` places it within ``max_subs`` of a
    human protein whose gene is transcribed above ``min_tpm`` in a tissue named by
    :data:`ESSENTIAL_TISSUES`. Reasons are
    ``[{"protein", "gene", "subs", "position", "tissue", "tpm"}]``.

    **This is a near-identity test, not a similarity test, and that distinction is the whole design.**
    The obvious alternative — score each register with :mod:`mhcmatch.mimicry` and flag the ones
    resembling a tolerance-side reference — was built and measured against this one on 1,000 viral
    epitopes (which cannot be self, so every firing is a false positive) and 1,000 thymic peptides
    from essential-tissue genes (``bench/results/vector_safety_screen.md``):

    ======================== ============ ============
    route                     false pos.   true pos.
    ======================== ============ ============
    mimicry, masked                0.693        0.944
    self origin, this one          0.020        0.940
    ======================== ============ ============

    **Equal sensitivity, 35× the false positives.** The reason is the one ``mimicry_collinear.md``
    already records: anchor-channel similarity to a presented reference is presentation, not
    recognition, so an anchor-masked match fires for every peptide sharing the allele's motif — the
    influenza epitope ``GILGFVFTL`` draws 14 essential-tissue hits. Nobody withdraws two-thirds of a
    cassette, so that route excludes nothing in practice and is not offered here.

    ``find_source`` separates instead: ``ESDPIVAQY`` resolves to ``sp|Q8WZ42|TITIN_HUMAN`` at 0
    substitutions, ``EVDPIGHLY`` to ``sp|P43357|MAGA3_HUMAN`` at 0 (and MAGE-A6 at 1), and
    ``GILGFVFTL`` to **nothing at all**. The 2 % that remains is not obviously noise — a viral 9-mer
    within one substitution of a human protein is what molecular mimicry means — and it is the floor
    on how specific this can get.

    **Two clauses, and the reason carries which one fired.** ``"target gene"`` — the unit's own
    ``gene`` is transcribed in an essential tissue, the MAGE-A12 case, and no register search is
    needed to see it. ``"unrelated self origin"`` — a register coincides within ``max_subs`` of a
    protein that is *not* the unit's own parent. Hits to the parent are dropped, because a long
    peptide is native context by design and tolerance already covers it; without that exclusion the
    screen rejects every unit of every cassette.

    **Clause 1 is skipped for a product whose sequence is not in the normal proteome, and that is a
    category error being corrected rather than a threshold being relaxed.** MAGE-A12 is a
    cancer-testis antigen: a shared, **unmutated** self protein, so its 0.33 TPM in brain caudate is
    the hazard exactly because the construct encodes a sequence brain tissue also presents. A somatic
    neoantigen is a different object — a missense, a frameshift, an inframe indel, a fusion junction
    all encode a sequence that **is absent from normal tissue by construction** — so the parent
    gene's expression is not that hazard. What is a hazard for it is clause 2, and clause 2 tests it
    for every kind, unchanged. Measured on a 37-donor cohort, clause 1 as an unconditional rule
    withdrew a candidate for the fact that its parent gene exists: **10 of 37 donors lost every unit
    they had**, and one lost 1,098 of 1,618 to clause 1 alone.

    ``novel_kinds`` is :data:`mhcmatch.predict.NOVEL_PRODUCTS` and is matched against
    :attr:`Unit.kind`. An ``isoform``, a ``cnv`` locus, a wild-type or overexpressed target is the
    MAGE-A12 case and **keeps** clause 1.

    **Clause 2 is asked only of the registers that carry novel sequence, and for the same reason.**
    A 27-mer unit is thirteen-fourteenths wild type by construction, and the unrestricted clause read
    that design as the hazard. Measured on **178 experimentally immunogenic somatic neoantigens**
    from ``isalgo/pmhc_data``, rebuilt as the 27-mer units they would enter a cassette as:
    **178 of 178 (100 %) trip clause 2**, at a median of **36** self registers each. 36 is
    ``12 + 10 + 8 + 6`` — exactly the count of 8/9/10/11-mer windows of a 27-mer that *cannot*
    contain a centred mutation — and the measured self fraction tracks that geometry at every length
    (L=8 60.02 % against 60.0 predicted, L=9 52.6 / 52.6, L=10 44.4 / 44.4, L=11 35.2 / 35.3;
    6,350 hits, **99.1 % of the geometric ceiling**). At the *minimal-epitope* level the clause is
    clean: **0 of 178** mutant epitopes are in the proteome and **178 of 178** wild types are. There
    were essentially no genuine coincidences to find — the veto was arithmetic, not evidence.

    So a window that does not contain novel sequence is **structurally exempt**: it is wild type, it
    was always going to be in the proteome, and no cassette avoids it short of not using long units.
    Which windows those are depends on the product, and
    :data:`mhcmatch.predict.TRACT_PRODUCTS` is the split: a ``frameshift`` or ``fusion`` is novel
    from :attr:`Unit.mutation_index` **to the end of the unit**, everything else in ``novel_kinds``
    at that one index. ``n_registers_spanning`` and ``n_hit_spanning`` ride on every clause-2 reason
    so the exemption is auditable rather than silent.

    The exemption is gated on the same ``novel_kinds`` as clause 1 — one list, two rules that cannot
    disagree. For an ``isoform``, a ``cnv`` or an unannotated unit **every** register is judged, as
    before, because for those the self-ness of the sequence *is* the finding. One case it reads
    generously: an in-frame fusion's downstream tract is genuine second-parent sequence, and it is
    exempted with the rest of the tract rather than charged to the second gene.

    **An unknown or empty kind keeps clause 1 — fail closed.** The screen may not exempt a unit
    because nobody annotated it. Every clause-1 reason therefore carries ``"kind"``, so a rejection
    can be read as *this is a shared self antigen* or as *nothing said what this was*. Note the one
    thing this cannot see: :attr:`Unit.kind` defaults to ``"missense"`` in the dataclass, so a
    :class:`Unit` constructed in Python with no ``kind`` is indistinguishable from one annotated as a
    somatic missense and is exempted. Annotating the unit is the caller's contract —
    :func:`units_from_context` and the CLI's unit table both fill it from the pipeline header.

    **``veto_tpm = 5.0`` separates a veto from a cost, and it is not the same line as ``min_tpm``.**
    ``min_tpm = 0.25`` stays what it always was: the *reporting* floor, set under MAGE-A12's 0.33 TPM
    so the fatal case is always visible. What 0.25 cannot also be is the *exclusion* line — at that
    level nearly every human gene is "detectable somewhere", which is what made the screen withdraw
    almost everything. ``veto_tpm`` is the conventional 5 TPM "is it expressed" cut, and with
    ``graded=True`` a finding below it is reported with ``"veto": False``: :func:`screen` keeps the
    unit, :func:`offtarget_cost` turns the finding into a per-unit cost, and
    :func:`mhcmatch.portfolio.compose` prices it against the response model instead of any one
    register vetoing a 27-mer. ``graded=False`` is the default and is the shipped veto behaviour.

    **``report_subs=1`` adds a third clause that reports and never withdraws.** The two deaths this
    screen is shaped around were both *near*-identity, not identity: titin's ``ESDPIVAQY`` differs
    from MAGE-A3's ``EVDPIGHLY`` at four positions, and MAGE-A12 is a different gene altogether. So
    the exact clause 2 cannot be the whole answer — but neither can a ``d=1`` veto, and the reason is
    measured rather than argued. On 178 validated immunogenic somatic neoantigens the exact clause
    withdraws **2 units (1.1%)**, while ``d=1`` to any different expressed gene reaches **125 units
    (70.2%)**: a veto there costs two thirds of every cassette to buy a hazard the exact clause has
    largely already taken. Clause 3 therefore emits ``"veto": False`` unconditionally — independently
    of ``graded`` — so :func:`screen` keeps the unit and ``notes`` carries the finding.

    Three filters keep that annotation readable, and the first carries most of it.
    ``report_min_length = 9`` excludes 8-mers, because at ``d = 1`` an 8-mer's ball is mostly chance:
    152 neighbours against 68,398,087 proteome windows in a space of 20**8 expects **0.41**
    coincidences per register, where a 9-mer's 171 neighbours in 20**9 expect **0.023** -- 18x fewer.
    Exact matching is unaffected and keeps its 8-mers, a d=0 8-mer expecting 0.0027 hits, which is
    why ``max_subs=0`` can scan a length ``report_subs=1`` must not. Then ``report_identity = 0.5``
    drops hits whose flanks are homologous to the unit's own context (:func:`flank_identity`),
    because a match to a paralog is descent rather than mimicry; and the off-target gene must clear
    ``min_tpm`` in an essential tissue, since a hazard needs something to be expressed.
    ``report_flank`` is how far either side the identity is read. Feed the survivors to
    :func:`presented` for the fourth and last filter. ``report_subs`` is 0 by default, which is
    exactly the two-clause screen as previously shipped.

    **``d=2`` is refused, not merely discouraged.** At radius 2 every expression floor from 0 TPM to
    100 TPM flags 178 of 178 units, with a median of 20 off-target genes each. The hazard genuinely
    does live out there — EPS8L2 at ``d=2``, titin at ``d=4`` — and no parameter reaches it without
    taking the entire cassette with it. That boundary is the finding, and it is why the screen stops
    at 1 and hands the rest to composition.

    **``max_subs=0`` — exact coincidence — because the decision is per unit while the search is per
    register, and that multiplies.** A 27-mer carries ~70 class-I registers and is withdrawn if any
    one of them fires, so a per-register false-positive rate that reads as small is not the rate a
    cassette experiences. Measured on six random 27-mers that carry no hazard, plus one burying the
    real titin epitope (``bench/results/vector_screen_radius.md``) — **units falsely withdrawn, out
    of the six**:

    ============ ======== ======== ==========
    ``max_subs`` 9-mers   9-11     8-11
    ============ ======== ======== ==========
    0            0        0        0
    1            1        1        **4**
    ============ ======== ======== ==========

    Radius 0 is clean at every length set. Radius 1 is not clean anywhere, and collapses once 8-mers
    enter: an 8-mer plus its 152 one-substitution neighbours is ~153 of 20\\ :sup:`8` sequences
    against the proteome's ~68 M windows, so a chance hit per register is expected and the ~20 8-mer
    registers in a unit make it near-certain. **Every setting still catches the titin unit**, so what
    radius 1 buys is nothing and what it costs is most of the cassette.

    **``min_tpm`` defaults to 0.25 because the two precedents disagree by two orders of magnitude and
    the lower one is what has to be caught.** Titin is 64.4 TPM in heart left ventricle and 351.4 in
    skeletal muscle, so any sane floor finds it. MAGE-A12 is **0.33 TPM** in brain caudate and 0.31 in
    putamen — expression that killed two patients, and that a conventional 5-TPM "is it expressed"
    cut would have waved through. The floor sits just under the fatal case, not at the conventional
    line. It still separates: MAGE-A3's own non-testis medians are 0.00.

    **What this does not catch, stated because a safety screen that oversells itself is worse than
    none.** It would not have caught the titin event *as it happened*. There the cassette contained
    MAGE-A3, whose profile is clean — 13.4 TPM in testis, 0.00 elsewhere — and the cross-reactive
    titin peptide was never in the construct. Four TCR-facing substitutions separate the two, so no
    distance threshold reaches it from the candidate, and the affinity-enhanced TCR that bridged them
    was the actual cause. What this catches is the **adjacent** and commoner failure: a register that
    *is* a self peptide from an essential-tissue gene, and a target gene like MAGE-A12 that is
    transcribed where it was assumed silent.

    ``symbols`` is ``{accession: gene}`` from
    :func:`mhcmatch.proteome.gene_symbols(path, key="accession")`; the search names proteins as
    ``sp|P43357|MAGA3_HUMAN`` and :func:`mhcmatch.expression.safety_profile` is keyed on ``MAGEA3``.
    It is required rather than defaulted because a missing map resolves nothing and so returns "no
    risk" for every peptide — the one wrong answer this must never give quietly.

    ``proteome`` needs :meth:`~mhcmatch.Proteome.find_sources`, the **batch** form, and at
    ``max_subs=0`` :meth:`~mhcmatch.Proteome.find_exact_sources` is used where the object has it: an
    exact question does not need the fuzzy index, which is a per-position Python loop costing ~12.6
    GB peak. The returned callable carries a ``prepare(registers)`` that :func:`screen` hands *every*
    register of *every* unit, so the whole candidate list resolves in one query rather than one per
    unit — screen everything in one process, never one unit per invocation.
    """
    from . import expression as EX

    if not symbols:
        raise ValueError("symbols must be a non-empty {accession: gene} map from "
                         "proteome.gene_symbols(path, key='accession'); without it every peptide "
                         "screens as safe")

    if report_subs > 1:
        raise ValueError(
            f"report_subs must be 0 (off) or 1, got {report_subs}. A radius of 2 is not a "
            "usable report tier: measured on 178 validated immunogenic neoantigens it flags "
            "178/178 at every expression floor, which annotates nothing.")

    novel = frozenset(str(k).strip().lower() for k in novel_kinds)
    tract = frozenset(str(k).strip().lower() for k in TRACT_PRODUCTS)
    prefixes = tuple(tissues)
    # {gene: [(tissue, tpm)]} above the reporting floor, computed once per gene. The filtered list
    # used to be rebuilt inside the loop -- once per unit for clause 1 and once per (register,
    # hit gene) for clause 2 -- so a cassette's worth of registers re-filtered the same handful of
    # genes thousands of times.
    ess: dict = {}
    # register -> [(gene, protein, n_subs, position)], one entry per DISTINCT source gene.
    cache: dict = {}

    def essential(gene):
        v = ess.get(gene)
        if v is None:
            # `top` must exceed the table's 123 distinct context names: at the default of 10 a
            # gene ranked by its highest tissues loses its essential-tissue rows before they are
            # ever tested, which is the second half of the same false negative as the naming fix.
            v = ess[gene] = [(t, x) for t, x in EX.safety_profile(gene, top=250)
                             if x >= min_tpm and t.startswith(prefixes)]
        return v

    # register -> [(gene, protein, position)] at exactly one substitution, own gene included and
    # nothing deduplicated: which of a gene's windows is the least homologous depends on the unit
    # asking, so the choice cannot be made here. ~1.06 hits per register, so the list stays short.
    near_cache: dict = {}

    def near(pep):
        """Every reference window exactly one substitution from ``pep``.

        The whole ``19 * len(pep)`` neighbourhood goes to :meth:`~mhcmatch.Proteome.find_exact_sources`
        as one batch, which answers it with a pair of ``searchsorted`` calls over the cached sorted
        window array. :meth:`~mhcmatch.Proteome.find_sources` with ``max_subs=1`` answers the same
        question by building the fuzzy index -- a Python loop over 68,389,335 positions -- and the
        two differ by more than three orders of magnitude on a cassette-sized batch.
        """
        v = near_cache.get(pep)
        if v is None:
            v = []
            for hs in proteome.find_exact_sources(_neighbours(pep)).values():
                for h in hs:
                    parts = h.protein.split("|")
                    gene = symbols.get(parts[1] if len(parts) >= 3 else h.protein)
                    if gene:
                        v.append((gene, h.protein, h.position))
            near_cache[pep] = v
        return v

    def novel_registers(unit, lengths):
        """The registers of ``unit`` that actually carry novel sequence, as a set.

        A point product (missense, an inframe indel, a stop/start loss) is novel at
        :attr:`Unit.mutation_index` and nowhere else, so a window is novel when it contains that
        index. A **tract** product (:data:`mhcmatch.predict.TRACT_PRODUCTS` -- frameshift, fusion)
        reads out of frame or across a junction from that index onward, so a window is novel when it
        reaches it at all.
        """
        mi, pep = int(getattr(unit, "mutation_index", 0)), unit.peptide
        wide = str(getattr(unit, "kind", "") or "").strip().lower() in tract
        # `i + L > mi` is "the window reaches the variant" -- the tract rule. A point product adds
        # "and does not start after it", which is containment.
        return {pep[i:i + L] for L in lengths for i in range(len(pep) - L + 1)
                if i + L > mi and (wide or i <= mi)}

    def prepare(registers):
        """Resolve every register's self origin in one query. Idempotent and additive."""
        todo = sorted({r for r in registers if r not in cache})
        if not todo:
            return
        exact = getattr(proteome, "find_exact_sources", None)
        hits = exact(todo) if (exact is not None and max_subs == 0) else \
            proteome.find_sources(todo, max_subs=max_subs)
        for pep in todo:
            seen, res = set(), []
            for h in hits.get(pep) or []:
                parts = h.protein.split("|")
                gene = symbols.get(parts[1] if len(parts) >= 3 else h.protein)
                if not gene or gene in seen:
                    continue
                seen.add(gene)
                res.append((gene, h.protein, h.n_subs, h.position))
            cache[pep] = res

    def risk(unit, registers):
        out = []
        kind = str(getattr(unit, "kind", "") or "").strip().lower()
        if kind not in novel:                               # clause 1: the target gene itself
            for tissue, tpm in essential(unit.gene):
                out.append({"clause": "target gene", "gene": unit.gene, "kind": kind,
                            "tissue": tissue, "tpm": tpm})
        # clause 2: an unrelated self origin, asked only of the registers that carry novel sequence.
        if kind in novel and registers and not any(r in unit.peptide for r in registers):
            # The spanning rule is defined against `unit.peptide`, so it is meaningless when the
            # registers are not that peptide's windows -- and it fails **silently**, by judging
            # nothing and reporting no hazard. That is the one answer a safety screen must not give
            # for a question it never asked. It happened: `bench/vector/safety_screen.py` handed a
            # bare placeholder unit the probe peptide and read 0 of 1000 on *both* arms, which reads
            # as a perfectly specific screen rather than a silent one.
            raise ValueError(
                f"none of the {len(registers)} registers passed is a window of unit.peptide, and "
                f"this unit's kind {kind!r} restricts clause 2 to the windows spanning "
                f"mutation_index {getattr(unit, 'mutation_index', '?')} -- so nothing would be "
                f"judged and the unit would screen as safe without being looked at. Pass the "
                f"unit's own registers, or set kind= outside NOVEL_PRODUCTS if the unit is a "
                f"shared, unmutated target.")
        judged = (novel_registers(unit, {len(r) for r in registers}) if kind in novel
                  else set(registers))
        registers = [r for r in registers if r in judged]
        prepare(registers)
        # The unit's own parent is native context, not a hazard: its flanks are self by construction
        # and its mutated register is one substitution from its own wild type.
        found = {pep: [h for h in cache.get(pep, ()) if h[0] != unit.gene] for pep in registers}
        counts = {"n_registers_spanning": len(registers),
                  "n_hit_spanning": sum(1 for v in found.values() if v)}
        for pep in registers:
            for gene, protein, subs, position in found[pep]:
                for tissue, tpm in essential(gene):
                    out.append({"clause": "unrelated self origin", "register": pep,
                                "protein": protein, "gene": gene, "subs": subs,
                                "position": position, "tissue": tissue, "tpm": tpm, **counts})
        if report_subs:
            # clause 3: a near-identical self origin in an unrelated gene. Never vetoes -- it is a
            # safety consideration attached to a kept unit, because a d=1 rule that withdrew would
            # cost two thirds of the candidates to buy a hazard the d=0 veto mostly already has.
            for pep in registers:
                if len(pep) < report_min_length:
                    continue
                at = unit.peptide.find(pep)
                best: dict = {}
                for gene, protein, position in near(pep):
                    if gene == unit.gene or not essential(gene):
                        continue
                    ident = flank_identity(unit.peptide, at, proteome.seqs[protein],
                                           position, len(pep), report_flank)
                    if ident < best.get(gene, (2.0,))[0]:
                        best[gene] = (ident, protein, position)
                for gene, (ident, protein, position) in sorted(best.items()):
                    if ident >= report_identity:
                        continue
                    tissue, tpm = max(essential(gene), key=lambda tv: tv[1])
                    out.append({"clause": "near-identical self origin", "register": pep,
                                "protein": protein, "gene": gene, "subs": report_subs,
                                "position": position, "identity": ident, "tissue": tissue,
                                "tpm": tpm, "veto": False,
                                "variant": proteome.seqs[protein][position:position + len(pep)],
                                **counts})
        if graded:
            for r in out:
                if r["tpm"] < veto_tpm:
                    r["veto"] = False
        return out

    risk.prepare = prepare
    return risk


def select(candidates, n0: float, cls: str | None = None, block=None) -> Selection:
    """Apply the per-allotype stopping rule (module docstring) to ranked ``candidates``.

    Candidates are grouped by :attr:`Unit.allele` and sorted by :attr:`Unit.p` descending within each
    group, then each group grows while ``p_next > S_a / (n0 + n_a)``. The first unit on an allotype is
    always taken: ``S_a = 0`` makes the threshold 0, and any candidate with ``p > 0`` clears it.

    ``n0`` is per-allotype capacity and must be positive. There is no default, deliberately — the
    literature does not fix it (see the module docstring), so a caller who has not decided what to
    assume has not finished designing the cassette.

    ``block`` chooses what the budget saturates against. The default is the allotype, which is the
    rule as shipped and as described in the module docstring. Passing a callable ``Unit -> hashable``
    blocks on something else, and the intended use is a *pair*: allotype together with the mechanism
    a unit was selected on, e.g.
    ``block=lambda u: (u.allele, corner[u.peptide])``. The arithmetic does not care what the key
    means; what it assumes is that two units sharing a key share a way of failing. Diversification
    across whatever that is falls out of the saturation, exactly as it does across allotypes — see
    :mod:`mhcmatch.portfolio` for the response model this is the greedy rule for, and for the
    measured intra-patient correlation that motivates blocking on more than the allotype.
    """
    if n0 <= 0:
        raise ValueError(f"n0 must be positive per-allotype capacity, got {n0}")
    key = block if block is not None else (lambda u: u.allele)
    pool = [c for c in candidates if cls is None or c.cls == cls]
    by = {}
    for c in pool:
        by.setdefault(key(c), []).append(c)

    kept, dropped, trace, keys = [], [], [], []
    for allele in sorted(by, key=str):
        ranked = sorted(by[allele], key=lambda u: (-u.p, u.gene, u.peptide))
        s, n = 0.0, 0
        stopped = False
        for c in ranked:
            threshold = s / (n0 + n)
            if stopped or c.p <= threshold:
                stopped = True                   # p is sorted descending: nothing later can clear it
                dropped.append(c)
                trace.append({"allele": allele, "gene": c.gene, "p": c.p,
                              "threshold": threshold, "kept": False})
                continue
            kept.append(c)
            keys.append(allele)
            trace.append({"allele": allele, "gene": c.gene, "p": c.p,
                          "threshold": threshold, "kept": True})
            s += c.p
            n += 1
    return Selection(units=kept, dropped=dropped, n0=float(n0), trace=trace, keys=keys)


def _path_cost(order_idx, cost) -> float:
    return sum(cost[order_idx[i]][order_idx[i + 1]] for i in range(len(order_idx) - 1))


def _greedy_2opt(cost, rounds: int = 4) -> list:
    """Cheapest-edge greedy open path, then 2-opt segment reversal until no improvement.

    Deterministic: no RNG anywhere, so the same units always give the same cassette. Simulated
    annealing (pVACvector's choice) explores more but cannot be reproduced from the inputs alone,
    which for a clinical artifact is the wrong trade.
    """
    n = len(cost)
    if n <= 2:
        return list(range(n))
    best = None
    for start in range(n):                       # every start, keep the cheapest -- n is <= tens
        unused, path = set(range(n)) - {start}, [start]
        while unused:
            nxt = min(unused, key=lambda j: (cost[path[-1]][j], j))
            path.append(nxt)
            unused.remove(nxt)
        if best is None or _path_cost(path, cost) < _path_cost(best, cost):
            best = path
    for _ in range(rounds):
        improved = False
        for i in range(n - 1):
            for j in range(i + 2, n):
                cand = best[:i + 1] + best[i + 1:j + 1][::-1] + best[j + 1:]
                if _path_cost(cand, cost) < _path_cost(best, cost) - 1e-12:
                    best, improved = cand, True
        if not improved:
            break
    return best


def order(units, binder, spacers=SPACERS, lengths=JUNCTION_LENGTHS,
          alleles=None, threshold: float | None = None,
          objective: str = "sum", binder_threshold: float | None = None) -> Cassette:
    """Choose a spacer and an ordering that minimise predicted junctional binding.

    **``objective`` matters and the two choices disagree, so it is explicit.**

    ``"sum"``   total of the strongest predicted binder at each junction. A junction is a hazard if
                *one* good binder forms there, which is pVACvector's logic (PMID 31907209). The
                junction count is ``n-1`` whatever the spacer, but a longer spacer creates more
                registers per junction and so a stochastically larger maximum — this objective
                therefore has a real bias toward the **shortest** spacer, up to and including none.
    ``"rate"``  predicted binders per register, needing ``binder_threshold``. Length-neutral, and it
                is the metric a junction sweep naturally reports.

    On one measured payload the two picked different spacers — ``"sum"`` chose no spacer where a
    rate sweep put ``AAA`` ahead of it — so a caller who has not chosen has not finished designing.

    Spacers are tried in ``spacers`` order and the **first** one whose worst junction falls at or
    below ``threshold`` wins; with ``threshold=None`` every spacer is tried and the one with the
    lowest total junction cost wins. Because :data:`SPACERS` leads with ``None``, a cassette that
    needs no spacer gets none — which is the right default, since every inserted residue is
    translated sequence that could itself form a binder.

    For each spacer the layout is an open path over a complete graph whose edge ``i -> j`` costs the
    strongest predicted binder spanning that junction, solved by :func:`_greedy_2opt`.

    ``binder(peptides, alleles) -> [float]``, higher meaning a stronger binder. Use
    :func:`store_binder` to build one from a :class:`~mhcmatch.store.Store`.
    """
    units = list(units)
    if len(units) < 2:
        seq = units[0].peptide if units else ""
        return Cassette(units=units, spacer=None, sequence=seq,
                        boundaries=[(0, len(seq))] if units else [], junctions=[], cost=0.0)

    best = None
    for sp in spacers:
        n = len(units)
        cost = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                pair = scan_junctions([units[i], units[j]], binder, sp, lengths, alleles)
                cost[i][j] = pair[0]["score"] if pair else float("-inf")
        idx = _greedy_2opt(cost)
        laid = [units[i] for i in idx]
        js = scan_junctions(laid, binder, sp, lengths, alleles, binder_threshold)
        finite = [j for j in js if j["score"] != float("-inf")]
        if objective == "rate":
            if binder_threshold is None:
                raise ValueError('objective="rate" needs binder_threshold')
            nw = sum(j["n_windows"] for j in js) or 1
            total = sum((j["n_over"] or 0) for j in js) / nw
        elif objective == "sum":
            total = sum(j["score"] for j in finite)
        else:
            raise ValueError(f'objective must be "sum" or "rate", got {objective!r}')
        worst = max((j["score"] for j in js), default=float("-inf"))
        cand = (total, worst, sp, laid, js)
        if best is None or total < best[0]:
            best = cand
        if threshold is not None and worst <= threshold:
            best = cand
            break

    total, _worst, sp, laid, js = best
    sequence, boundaries, at = "", [], 0
    for k, u in enumerate(laid):
        if k:
            sequence += sp or ""
            at += len(sp or "")
        boundaries.append((at, at + len(u.peptide)))
        sequence += u.peptide
        at += len(u.peptide)
    return Cassette(units=laid, spacer=sp, sequence=sequence, boundaries=boundaries,
                    junctions=js, cost=float(total))


#: Codons synonymous with ``TTT`` (phenylalanine). ``TTC`` is the only alternative, which is what
#: makes :func:`deslip` a one-line fix rather than a codon-optimisation problem.
_PHE = ("TTT", "TTC")


def slippery_sites(cds: str) -> list:
    """Codon positions where N1-methylpseudouridine drives **+1 ribosomal frameshifting**.

    Returns ``[{codon_index, nt_offset, codon, next_codon}, ...]``.

    m1Ψ is not translationally neutral. Mulroney et al. (*Nature* 2024, PMID 38057663) measured **+1
    frameshifting at ~8% of the in-frame product** in m1Ψ mRNA, localised to a slippery motif —
    ``m1Ψ m1Ψ m1Ψ X`` with ``X`` = m1Ψ or C at the **first position of the following codon**, i.e. a
    ``TTT`` codon followed by a codon starting ``T`` or ``C``. **Six such sites sit in the BNT162b2
    spike coding sequence**, and BNT162b2-vaccinated humans mounted a significantly higher IFN-γ
    response against the +1 frameshifted product than controls.

    **This matters far more for a designed polyepitope than for a natural ORF**, for two reasons.
    A concatemer has many more codon-boundary junctions per kilobase, and the residues at those
    junctions are the designer's choice — glycine/serine linkers are encoded by ``GGN``/``AGY``/``TCN``,
    which is exactly how U-runs end up at seams. And the consequence is worse: a frameshift inside a
    polyepitope does not merely lose protein, it translates an entire downstream out-of-frame
    cassette that is itself presented, so the construct delivers a second, unintended and unscreened
    antigen payload.

    Scanning is therefore **mandatory for any m1Ψ construct** and pointless for an unmodified-uridine
    one — BioNTech's cancer platform deliberately uses unmodified uridine, so which applies depends
    on the platform, not the sequence.

    Only the codon-aligned motif as published is reported. Whether non-codon-aligned U-runs also
    induce slippage was not characterised in that work, so it is not guessed at here.
    """
    s = cds.strip().upper().replace("U", "T")
    out = []
    for i in range(0, len(s) - 5, 3):
        codon, nxt = s[i:i + 3], s[i + 3:i + 6]
        if codon == "TTT" and nxt[:1] in ("T", "C"):
            out.append({"codon_index": i // 3, "nt_offset": i,
                        "codon": codon, "next_codon": nxt})
    return out


def deslip(cds: str) -> tuple:
    """``(cds, n_fixed)`` — remove every :func:`slippery_sites` motif synonymously.

    ``TTT`` and ``TTC`` both encode phenylalanine, so rewriting the *upstream* codon breaks the U-run
    without touching the protein and without disturbing the downstream codon's own optimisation.
    This is the fix Mulroney et al. validated: single ``U*187C`` / ``U*208C`` substitutions strongly
    reduced frameshifting and the double mutant produced none detectable.

    Idempotent — a second call finds nothing to do.
    """
    s = list(cds.strip().upper().replace("U", "T"))
    n = 0
    for site in slippery_sites("".join(s)):
        i = site["nt_offset"]
        s[i:i + 3] = list("TTC")
        n += 1
    return "".join(s), n


#: Homo sapiens codon usage as ``{codon: (amino acid, occurrences per thousand codons)}``.
#: Kazusa Codon Usage Database (https://www.kazusa.or.jp/codon/), *Homo sapiens* [gbpri] --
#: 93,487 CDSs, 40,662,582 codons; retrieved 2026-08-18.
#:
#: Stored as the **measured frequencies**, not as a pre-reduced "best codon per residue" map, so
#: :func:`back_translate`'s choice is derived from data under a stated rule and can be re-derived
#: under a different one. Substituting a host's own table is then a one-argument change.
CODON_USAGE_HUMAN: dict = {
    "TTT": ("F", 17.6), "TTC": ("F", 20.3), "TTA": ("L", 7.7),  "TTG": ("L", 12.9),
    "CTT": ("L", 13.2), "CTC": ("L", 19.6), "CTA": ("L", 7.2),  "CTG": ("L", 39.6),
    "ATT": ("I", 16.0), "ATC": ("I", 20.8), "ATA": ("I", 7.5),  "ATG": ("M", 22.0),
    "GTT": ("V", 11.0), "GTC": ("V", 14.5), "GTA": ("V", 7.1),  "GTG": ("V", 28.1),
    "TCT": ("S", 15.2), "TCC": ("S", 17.7), "TCA": ("S", 12.2), "TCG": ("S", 4.4),
    "AGT": ("S", 12.1), "AGC": ("S", 19.5),
    "CCT": ("P", 17.5), "CCC": ("P", 19.8), "CCA": ("P", 16.9), "CCG": ("P", 6.9),
    "ACT": ("T", 13.1), "ACC": ("T", 18.9), "ACA": ("T", 15.1), "ACG": ("T", 6.1),
    "GCT": ("A", 18.4), "GCC": ("A", 27.7), "GCA": ("A", 15.8), "GCG": ("A", 7.4),
    "TAT": ("Y", 12.2), "TAC": ("Y", 15.3), "TAA": ("*", 1.0),  "TAG": ("*", 0.8),
    "CAT": ("H", 10.9), "CAC": ("H", 15.1), "CAA": ("Q", 12.3), "CAG": ("Q", 34.2),
    "AAT": ("N", 17.0), "AAC": ("N", 19.1), "AAA": ("K", 24.4), "AAG": ("K", 31.9),
    "GAT": ("D", 21.8), "GAC": ("D", 25.1), "GAA": ("E", 29.0), "GAG": ("E", 39.6),
    "TGT": ("C", 10.6), "TGC": ("C", 12.6), "TGA": ("*", 1.6),  "TGG": ("W", 13.2),
    "CGT": ("R", 4.5),  "CGC": ("R", 10.4), "CGA": ("R", 6.2),  "CGG": ("R", 11.4),
    "AGA": ("R", 12.2), "AGG": ("R", 12.0),
    "GGT": ("G", 10.8), "GGC": ("G", 22.2), "GGA": ("G", 16.5), "GGG": ("G", 16.5),
}

#: Run of one nucleotide past which :func:`back_translate` reaches for a rarer synonymous codon.
#: Homopolymers are a **synthesis** constraint, not a translation one: vendors reject or
#: mis-assemble long single-base runs, and a spacered concatemer manufactures them directly
#: (``AAA`` and ``GPGPG`` are both in :data:`SPACERS`). Four is the conventional screening floor.
#:
#: It is a **target, not a guarantee**, because the choice is greedy and per-codon. Measured over
#: 5,000 random 20-60mers under the default table: longest run 6, and 84% of sequences at or below
#: 4; the same peptides back-translated by most-frequent-codon alone reach 13, with 6% above 6.
MAX_HOMOPOLYMER: int = 4


def _synonyms(usage: dict) -> dict:
    """``{amino acid: (codon, ...)}`` ordered by descending usage, then codon for determinism."""
    out: dict = {}
    for codon, (aa, freq) in usage.items():
        out.setdefault(aa, []).append((freq, codon))
    return {aa: tuple(c for _, c in sorted(v, key=lambda t: (-t[0], t[1]))) for aa, v in out.items()}


def _longest_run(seq: str) -> int:
    """Longest run of one repeated nucleotide anywhere in ``seq``."""
    best = run = 0
    last = ""
    for ch in seq:
        run = run + 1 if ch == last else 1
        last = ch
        best = max(best, run)
    return best


def translate(cds: str) -> str:
    """Amino-acid sequence of a coding sequence, ``*`` for a stop.

    Exists so "synonymous" is checkable rather than asserted -- :func:`deslip` and
    :func:`back_translate` both claim it, and a caller who supplies their own codon table needs the
    same check. ``U`` is read as ``T``; a trailing partial codon is ignored.
    """
    s = cds.strip().upper().replace("U", "T")
    return "".join(CODON_USAGE_HUMAN.get(s[i:i + 3], ("X", 0.0))[0] for i in range(0, len(s) - 2, 3))


def back_translate(peptide: str, usage: dict = None, *, avoid_slip: bool = True,
                   max_run: int = MAX_HOMOPOLYMER) -> str:
    """Coding sequence for ``peptide`` -- the cassette's nucleotide half.

    **Highest-usage synonymous codon per residue**, backing off to the next one whenever the first
    would extend a single-nucleotide run past ``max_run``, then :func:`deslip` to remove the m1-psi
    +1-frameshift motif. Deterministic: the same peptide and table always give the same CDS.

    The backoff is greedy, so ``max_run`` is a target rather than a bound -- see
    :data:`MAX_HOMOPOLYMER` for what it is measured to be worth.

    What this is *not* is a codon optimiser. It fixes the two things that make a **polyepitope**
    construct fail where a natural ORF would not -- the frameshift motif (:func:`slippery_sites`,
    which a concatemer hits far more often because the designer chooses the seam residues) and
    synthesis-hostile homopolymers (which spacers like ``AAA`` manufacture directly). It does not
    touch GC content, secondary structure, splice sites or CpG, and a manufacturer's own optimiser
    should be preferred where one is available; this exists so a cassette ships with a usable CDS
    rather than none.

    Emits the epitope cassette only -- no start codon, no stop, no leader, no trafficking domain --
    matching :attr:`Cassette.sequence`, because those flanks are the vector's, not the payload's.

    >>> cds = back_translate("SIINFEKL")
    >>> translate(cds)
    'SIINFEKL'
    >>> slippery_sites(cds)
    []
    """
    table = usage or CODON_USAGE_HUMAN
    syn = _synonyms(table)
    out: list = []
    seq = ""
    for aa in peptide.strip().upper():
        options = syn.get(aa)
        if not options:
            raise ValueError(f"no codon for residue {aa!r} in the supplied usage table")
        # The run a codon creates can sit *inside* it, so the whole junction window is checked and
        # not just the trailing run -- ...AAA + AAG ends in G but carries AAAAA across the seam.
        # Two tiers: the most-used codon that stays within ``max_run``, and failing that the one
        # that makes the shortest run at all. The second tier is not hypothetical -- proline's four
        # codons all begin ``CC``, so consecutive prolines cannot be brought below a 5-run by any
        # synonymous choice, and the rule degrades to "as short as this residue allows".
        tail = seq[-max_run:]
        runs = [(_longest_run(tail + c), i, c) for i, c in enumerate(options)]
        within = [t for t in runs if t[0] <= max_run]
        pick = min(within or runs, key=lambda t: t[1] if within else (t[0], t[1]))[2]
        out.append(pick)
        seq += pick
    cds = "".join(out)
    return deslip(cds)[0] if avoid_slip else cds


def store_binder(store, alleles, cls: str = "mhc1"):
    """A ``binder`` callable over a :class:`~mhcmatch.store.Store`: ``-log10(%rank)`` of the best
    allele, so higher means a stronger predicted binder.

    Kept as a thin adapter, and out of :func:`order`, so the layout logic stays testable with no
    panel, no download and no allele list.
    """
    import math

    def binder(peptides, _alleles=None):
        out = []
        for p in peptides:
            r = store.restriction(p, cls=cls, alleles=list(_alleles or alleles), calibrated=True)
            # Restriction.rank is the per-allele %rank and is None unless calibrated=True; a peptide
            # with no restriction at all is a non-binder, i.e. %rank 100.
            pct = min((x.rank for x in r if x.rank is not None), default=100.0)
            out.append(-math.log10(max(pct, 1e-4)))
        return out

    return binder


def rebuild(cassette: Cassette, **kw) -> Cassette:
    """Re-lay an existing cassette's units under different settings, keeping the units fixed.

    The point of comparison is the incumbent: re-order and re-space the *same* payload to separate
    "these units were the wrong choice" from "these units were laid out badly".
    """
    return order(cassette.units, **kw)


def from_sequence(sequence: str, spacer: str, lengths=JUNCTION_LENGTHS) -> list:
    """Split an existing cassette on a known ``spacer`` into pseudo-units, for auditing.

    Only usable when the spacer is unambiguous — a cassette joined by alternating tokens has to be
    split by the caller, who knows the grammar. ``mutation_index`` is unknown from sequence alone and
    is set to the window centre; ``p`` is 0.0. These units are for junction scanning, not selection.
    """
    parts = [s for s in sequence.split(spacer) if s]
    return [Unit(peptide=s, mutation_index=len(s) // 2, gene=f"seg{k}", allele="", p=0.0)
            for k, s in enumerate(parts)]


# ---------------------------------------------------------------------------- the cassette map

#: Class-II ligand lengths scanned when mapping a cassette. Wider than
#: :data:`MHC2_JUNCTION_LENGTHS`, which exists to score *junctions* and only needs the shortest
#: windows a core can be read from: a map is annotating what a cassette actually presents, and a
#: class-II ligand runs to 25.
MHC2_MAP_LENGTHS: tuple = tuple(range(12, 21))


@dataclass(frozen=True)
class Feature:
    """One annotated span of an assembled cassette, in **1-based inclusive** amino-acid coordinates.

    ``kind`` is ``unit`` (a vaccine unit), ``linker`` (the spacer between two of them) or
    ``epitope`` (a predicted binder). Units and linkers tile the cassette exactly; epitopes overlay
    it and may span a junction, which is the case ``unit = 0`` marks.
    """

    id: str
    kind: str
    start: int
    end: int
    seq: str
    cls: str = ""
    allele: str = ""
    rank: float = float("nan")
    unit: int = 0
    gene: str = ""
    core_start: int = 0
    core_end: int = 0
    #: The binding core itself (:func:`mhcmatch.store.binding_core`), 9 residues, both classes.
    #: ``core_start``/``core_end`` stay class-II-only because a class-I core is not a
    #: contiguous span -- it drops the bulge -- so it has no cassette coordinates to give.
    core: str = ""
    overlaps: tuple = ()

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def store_ranker(store, alleles, cls: str = "mhc1", calibrated: bool = True):
    """A ``ranker`` callable over a :class:`~mhcmatch.store.Store`: ``[(allele, %rank), ...]`` per
    peptide, **one entry per allele that presents it**.

    Distinct from :func:`store_binder`, which collapses to the best allele because a *layout* cost
    only needs to know whether some binder forms. A map needs the allele: a heterozygote presents
    the same peptide on two molecules, and that is two facts about the cassette rather than one.
    """
    def ranker(peptides):
        out = []
        for p in peptides:
            r = store.restriction(p, cls=cls, alleles=list(alleles), calibrated=True)
            out.append([(x.allele, float(x.rank)) for x in r if x.rank is not None])
        return out
    return ranker


def _windows(seq: str, lengths) -> list:
    return [(i, i + L, seq[i:i + L]) for L in lengths for i in range(len(seq) - L + 1)]


def epitope_map(cassette: Cassette, ranker1=None, ranker2=None, threshold: float = 2.0,
                lengths1=JUNCTION_LENGTHS, lengths2=MHC2_MAP_LENGTHS) -> list:
    """Annotate an assembled cassette: units, linkers, predicted epitopes, and **which class-I and
    class-II epitopes overlap each other**.

    **Why the overlap is the point and not a decoration.** A cassette that carries a CD8 epitope and
    borrows its CD4 help from an unrelated universal helper (PADRE, HBVcore) raises no T-cell
    response against the tumour antigen on the class-II side. Kissick *et al.* built one 27-mer
    around the HLA-A\\*02:01 SIM2\\ :sub:`237-245` epitope so that a class-II epitope from the **same**
    protein overlapped it, and it replaced the exogenous HBVcore helper outright: the long peptide
    alone raised both the CD8 IFN-γ recall response to the 9-mer and a CD4 IL-2 response to
    SIM2\\ :sub:`240-254`, and 137 class-II binders were predicted across DR/DP/DQ from that one
    27-mer (*PLoS One* 2014;9(4):e93231, PMID 24690990, doi:10.1371/journal.pone.0093231). A unit
    whose class-I epitope has **no** overlapping class-II epitope is the configuration that needed
    the borrowed helper, and this map is what says which units those are.

    ``ranker1`` / ``ranker2`` are ``ranker(peptides) -> [[(allele, %rank), ...], ...]``, per class;
    :func:`store_ranker` builds one from a :class:`~mhcmatch.store.Store`. Either may be ``None``,
    which simply omits that class. Injected for the same reason :func:`order`'s ``binder`` is — the
    whole map is testable with no panel, no download and no calibration.

    Every ``(peptide, allele)`` at or below ``threshold`` %rank is its own :class:`Feature`, so a
    peptide presented by two of the patient's alleles appears **twice**. That is deliberate: at a
    heterozygous locus the two molecules are two independent presentation events, they are what the
    per-allotype capacity in :func:`select` is spent on, and collapsing them would under-count the
    cassette's coverage of exactly the patients it was personalised for.

    Coordinates are 1-based inclusive over :attr:`Cassette.sequence`, which is the epitope cassette
    only — no start codon, no leader, no tag. An mRNA construct that adds those must offset.
    """
    seq = cassette.sequence
    feats: list = []

    # Units and linkers tile the cassette; `boundaries` is 0-based half-open.
    prev_end = 0
    for i, ((lo, hi), u) in enumerate(zip(cassette.boundaries, cassette.units), 1):
        if lo > prev_end:
            feats.append(Feature(id=f"l{i - 1}", kind="linker", start=prev_end + 1, end=lo,
                                 seq=seq[prev_end:lo]))
        feats.append(Feature(id=f"u{i}", kind="unit", start=lo + 1, end=hi, seq=seq[lo:hi],
                             unit=i, gene=u.gene, allele=u.allele or ""))
        prev_end = hi

    def in_unit(lo, hi):
        """Which unit fully contains [lo, hi); 0 when it spans a junction or a linker."""
        for i, (a, b) in enumerate(cassette.boundaries, 1):
            if a <= lo and hi <= b:
                return i
        return 0

    n = 0
    for cls, ranker, lengths in (("mhc1", ranker1, lengths1), ("mhc2", ranker2, lengths2)):
        if ranker is None:
            continue
        wins = _windows(seq, lengths)
        if not wins:
            continue
        for (lo, hi, pep), hits in zip(wins, ranker([w[2] for w in wins])):
            for allele, rank in hits:
                if rank is None or rank > threshold:
                    continue
                n += 1
                k = in_unit(lo, hi)
                from .store import binding_core
                core_seq, off = binding_core(pep, cls)
                span = (lo + off + 1, lo + off + 9) if (cls == "mhc2" and core_seq) else (0, 0)
                feats.append(Feature(id=f"e{n}", kind="epitope", start=lo + 1, end=hi, seq=pep,
                                     cls=cls, allele=allele, rank=float(rank), unit=k,
                                     gene=cassette.units[k - 1].gene if k else "",
                                     core_start=span[0], core_end=span[1], core=core_seq))

    # Cross-class overlap, computed once over the finished list so both directions agree.
    e1 = [f for f in feats if f.kind == "epitope" and f.cls == "mhc1"]
    e2 = [f for f in feats if f.kind == "epitope" and f.cls == "mhc2"]
    over: dict = {}
    for a in e1:
        for b in e2:
            if a.start <= b.end and b.start <= a.end:
                over.setdefault(a.id, []).append(b.id)
                over.setdefault(b.id, []).append(a.id)
    return [f if f.id not in over else replace(f, overlaps=tuple(over[f.id])) for f in feats]


#: Column order of the cassette map, one source of truth for the TSV and the JSON.
MAP_COLUMNS: tuple = ("id", "kind", "start", "end", "length", "seq", "cls", "allele", "rank",
                      "unit", "gene", "core", "core_start", "core_end", "overlaps", "n_overlaps")


def map_rows(features) -> list:
    """The map as plain dicts in :data:`MAP_COLUMNS` order — the payload both writers share."""
    out = []
    for f in features:
        out.append({"id": f.id, "kind": f.kind, "start": f.start, "end": f.end, "length": f.length,
                    "seq": f.seq, "cls": f.cls, "allele": f.allele,
                    "rank": None if f.rank != f.rank else round(f.rank, 4),
                    "unit": f.unit, "gene": f.gene, "core": f.core,
                    "core_start": f.core_start or None, "core_end": f.core_end or None,
                    "overlaps": list(f.overlaps), "n_overlaps": len(f.overlaps)})
    return out


def map_summary(cassette: Cassette, features) -> dict:
    """Per-unit coverage: how many class-I and class-II epitopes, on how many allotypes, and
    **whether the unit's class-I epitopes have class-II help from within the same unit**.

    The last column is the one a reviewer reads first. A unit with class-I epitopes and no
    overlapping class-II epitope is the configuration that needed a borrowed universal helper.
    """
    eps = [f for f in features if f.kind == "epitope"]
    units = []
    for i, u in enumerate(cassette.units, 1):
        mine = [f for f in eps if f.unit == i]
        c1 = [f for f in mine if f.cls == "mhc1"]
        c2 = [f for f in mine if f.cls == "mhc2"]
        helped = [f for f in c1 if f.overlaps]
        units.append({
            "unit": i, "gene": u.gene, "allele": u.allele or "", "peptide": u.peptide,
            "start": cassette.boundaries[i - 1][0] + 1, "end": cassette.boundaries[i - 1][1],
            "n_mhc1": len(c1), "n_mhc2": len(c2),
            "alleles_mhc1": sorted({f.allele for f in c1}),
            "alleles_mhc2": sorted({f.allele for f in c2}),
            "n_mhc1_with_mhc2_overlap": len(helped),
            "self_help": bool(helped),
        })
    spanning = [f for f in eps if f.unit == 0]
    return {
        "n_units": len(cassette.units), "spacer": cassette.spacer,
        "length_aa": len(cassette.sequence),
        "n_mhc1": sum(f.cls == "mhc1" for f in eps),
        "n_mhc2": sum(f.cls == "mhc2" for f in eps),
        "n_junction_spanning": len(spanning),
        "n_units_with_self_help": sum(u["self_help"] for u in units),
        "units": units,
    }


def write_map(cassette: Cassette, features, tsv_path: str | None = None,
              json_path: str | None = None) -> dict:
    """Write the cassette map as TSV and/or JSON, and return the summary dict.

    The TSV is the flat table — one row per feature, one value per cell, so it sorts and joins. The
    JSON carries the same rows **plus** the per-unit summary and the cassette sequence, which is what
    a viewer needs to draw the thing without recomputing anything.
    """
    import csv
    import json as _json

    rows = map_rows(features)
    summary = map_summary(cassette, features)
    if tsv_path:
        with open(tsv_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(MAP_COLUMNS), delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow({**r, "overlaps": ",".join(r["overlaps"])})
    if json_path:
        with open(json_path, "w") as fh:
            _json.dump({"sequence": cassette.sequence, "spacer": cassette.spacer,
                        "summary": summary, "features": rows}, fh, indent=1)
    return summary
