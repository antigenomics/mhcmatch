"""Assembling a polyepitope vaccine cassette: how many units, in what order, joined by what.

Candidate *selection* — which mutations are worth targeting — is :mod:`mhcmatch.rank` and the
immunogenicity stack behind it. This module is the step after: given ranked candidates with
calibrated probabilities, decide **how many to carry, how to lay them out, and what to put between
them**. Those are three separate questions with three different literatures, and none of them is
answered by the candidate score.

Selection and assembly are kept apart because the assembly answer depends on the *set*, not the
candidate: whether to carry a 12th epitope depends on what the first eleven already cover, and the
cost of a junction depends on which two units sit either side of it.

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

``AAY``    ends in tyrosine. ERAP1 prefers hydrophobic C-termini and has low affinity for charged
           ones (Chang et al., *PNAS* 2005, PMID 16286653), so a terminal Y genuinely aids
           processing — while also supplying the C-terminal anchor for A\\*01:01, A\\*29:02 and
           B\\*35:01. It is a trade-off, not a mistake, and which way it falls is donor-specific.
``KK``     leaves charged residues at the boundary: the mirror image, poor for ERAP1.
G/P-rich   glycine and proline are disfavoured at MHC-I anchor positions and abundant in the
           C-terminal regions from which ligands are cleaved (Martín-Galiano & López, *PLoS One*
           2019, PMID 30645615), so these sit in the permissive zone.

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
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

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

    @property
    def expected_yield(self) -> float:
        """``sum_a n0 * S_a / (n0 + n_a)`` — expected responses under the saturation model."""
        by = {}
        for u in self.units:
            by.setdefault(u.allele, []).append(u.p)
        return sum(self.n0 * sum(ps) / (self.n0 + len(ps)) for ps in by.values())

    def per_allele(self) -> dict:
        """``{allele: (n_units, summed p, saturated yield)}`` — where the budget actually went."""
        by = {}
        for u in self.units:
            by.setdefault(u.allele, []).append(u.p)
        return {a: (len(ps), sum(ps), self.n0 * sum(ps) / (self.n0 + len(ps)))
                for a, ps in sorted(by.items())}


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
                   lengths=JUNCTION_LENGTHS, alleles=None) -> list:
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
        out.append({"left": k, "right": k + 1, "score": float(scores[b]),
                    "peptide": peps[b], "offset": wins[b][1], "n_windows": len(peps)})
    return out


def select(candidates, n0: float, cls: str | None = None) -> Selection:
    """Apply the per-allotype stopping rule (module docstring) to ranked ``candidates``.

    Candidates are grouped by :attr:`Unit.allele` and sorted by :attr:`Unit.p` descending within each
    group, then each group grows while ``p_next > S_a / (n0 + n_a)``. The first unit on an allotype is
    always taken: ``S_a = 0`` makes the threshold 0, and any candidate with ``p > 0`` clears it.

    ``n0`` is per-allotype capacity and must be positive. There is no default, deliberately — the
    literature does not fix it (see the module docstring), so a caller who has not decided what to
    assume has not finished designing the cassette.
    """
    if n0 <= 0:
        raise ValueError(f"n0 must be positive per-allotype capacity, got {n0}")
    pool = [c for c in candidates if cls is None or c.cls == cls]
    by = {}
    for c in pool:
        by.setdefault(c.allele, []).append(c)

    kept, dropped, trace = [], [], []
    for allele in sorted(by):
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
            trace.append({"allele": allele, "gene": c.gene, "p": c.p,
                          "threshold": threshold, "kept": True})
            s += c.p
            n += 1
    return Selection(units=kept, dropped=dropped, n0=float(n0), trace=trace)


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
          alleles=None, threshold: float | None = None) -> Cassette:
    """Choose a spacer and an ordering that minimise predicted junctional binding.

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
        js = scan_junctions(laid, binder, sp, lengths, alleles)
        total = sum(j["score"] for j in js if j["score"] != float("-inf"))
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
            pct = min((x.percent_rank for x in r), default=100.0) if r else 100.0
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
