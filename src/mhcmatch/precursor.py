"""T-cell precursor frequency for an epitope: how much repertoire mass can see it.

The estimand is

    F(e) = sum over the cognate set C_e of pi(tau)

— the probability that a random naive-repertoire junction recognises epitope ``e``. This is the
continuous quantity behind "immunogenic": recognition is a spectrum, and `F(e)` is where on it a
given epitope sits.

Five estimators of the same `F(e)`, in increasing order of model commitment:

``observed_mass``
    Sum of `Pgen` over the junctions actually recorded for the epitope. A **strict lower bound**,
    and a biased one: VDJdb samples cognate TCRs **size-biased by Pgen** (a TCR enters the record
    roughly in proportion to its repertoire frequency), so the observed members are systematically
    the high-`Pgen` ones and the deficit does not shrink with more studies at the same depth.

``coverage_corrected_mass``
    ``observed_mass`` with the size-biased deficit put back. A capture curve increasing in `Pgen` is
    fitted to the per-junction donor/study multiplicities, then the mass is Horvitz–Thompson
    reweighted by each member's inclusion probability. Unlike the bound it is *not* monotone in
    depth by construction, and it degenerates loudly rather than silently — see the function's own
    docstring for what it assumes and where it breaks.

``ball_mass``
    Mass of the **union** of Hamming-`r` balls around the observed junctions. Cognate TCRs are
    near-duplicates by construction, so their balls overlap and the *sum* of per-sequence ball
    masses double-counts; the union is the correct object. The returned ``overlap`` quantifies that
    double-counting directly, which is worth reporting rather than hiding — it is a measurement of
    how tight the specificity group is.

``shell_profile``
    ``ball_mass`` resolved by exact edit distance, so the empirically measured cognacy-retention
    profile `alpha_r` can be applied per shell instead of assuming every ball member is still
    cognate. `F ≈ sum_r alpha_r * mass(shell r)`.

``motif_mass``
    `Pgen` of a degenerate motif (a VDJdb cluster PWM is V/J/length-pinned, hence exactly a
    per-position residue set — see :func:`load_cluster_motifs`). This computes the mass of a **set**
    directly, so unlike the others it does not suffer the observed-sample coverage bias at all.

Which one to use
----------------

============================================  =======================================
question                                      estimator
============================================  =======================================
"what can I defend without any assumption?"   ``observed_mass`` (report it as a bound)
"how much did the sampling miss?"             ``coverage_corrected_mass`` (needs >= 2
                                              capture units and some recaptures)
"how tight is this specificity group?"        ``ball_mass`` -> the ``overlap`` field
"best point estimate from observed TCRs"      ``shell_profile`` -> ``retained``
"an estimate with no sampling bias at all"    ``motif_mass`` (needs a cluster PWM)
"how much mass is still missing?"             ``cross_check`` -> ``missing_fraction``
============================================  =======================================

:func:`cross_check` is the scientifically load-bearing one. ``motif_mass`` and the observed-sample
estimators measure the same `F(e)` by independent routes with *different* biases, so **their
disagreement is an estimate of the missing mass**.

Requires the optional ``vdjtools`` dependency (the recombination model). Nothing here reimplements
`Pgen` — the DP, the closed Hamming-1 ball and the degenerate/masked DP all live in vdjtools, and
the neighbourhood enumeration lives in seqtree.
"""

from __future__ import annotations

import csv
import gzip
import math
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "load_model", "check_junctions", "pgen",
    "observed_mass", "coverage_corrected_mass",
    "ball_mass", "shell_profile",
    "ClusterMotif", "load_cluster_motifs", "motif_mass",
    "cross_check",
    "ALPHA_PER_EDIT", "MOTIF_FREQ_THRESHOLD", "MAX_BALL_MEMBERS",
]

_JUNCTION_START = "C"
_JUNCTION_END = ("F", "W")

#: Cognacy retention per unit of edit distance -- the probability that a junction one substitution
#: away from a cognate TCR is itself cognate. Measured by Mayer & Callan, *PNAS* 2023;120:e2213264120
#: (PMID 36649423) from near-coincidence statistics within epitope-specific repertoires: binding
#: probability falls roughly **ten-fold per Levenshtein unit**, consistently across disparate
#: experiments. It is a parameter of :func:`shell_profile`, not a constant of the library -- pass
#: your own if you have measured one for your data.
ALPHA_PER_EDIT = 0.1

#: Default per-position residue-frequency cut for turning a VDJdb cluster PWM into an allowed-set
#: motif. Calibrated on the 2026-06 release against ``cluster_members.txt``: the largest round
#: threshold at which the motifs still match >= 95% of their own clusters' member junctions
#: (measured 0.10 -> 99.2% member recall, 0.15 -> 96.6%, 0.20 -> 92.1%). Raise it for a tighter,
#: lower-mass motif; lower it towards 0 for the union of everything observed in the cluster.
MOTIF_FREQ_THRESHOLD = 0.15

#: Ceiling on the number of enumerated neighbourhood members, for :func:`shell_profile`. The
#: enumeration is materialised as Python strings: 300 measured junctions at ``r=2`` produce ~9.9M
#: sequences and cost ~1.8 GB, i.e. ~190 bytes each. The default keeps a single call under ~0.4 GB.
#: ``r=1`` is ~19L per junction and never comes close; ``r=2`` is ~180 L^2/2 and usually does.
MAX_BALL_MEMBERS = 2_000_000


def _native():
    """The vdjtools Pgen bindings, or a clear error naming the extra."""
    try:
        from vdjtools.model import native
    except ImportError as e:                     # pragma: no cover - depends on the environment
        raise ImportError(
            "mhcmatch.precursor needs the optional 'vdjtools' dependency for the recombination "
            "model (Pgen). Install it with: pip install 'mhcmatch[precursor]'"
        ) from e
    return native


def load_model(locus: str = "TRB", source: str = "olga", organism: str = "human"):
    """Bundled recombination model.

    ``source="olga"`` is the default for null work. Do **not** use ``"learned"``: those models were
    EM-fit on ~2k clonotypes without a gene-usage pseudocount, so 68 of 89 bundled TRB V alleles
    have `P(V) = 0` and any junction using them scores 0. Mouse is available only under
    ``source="arda"``, and only for TRA/TRB.
    """
    from vdjtools.model import load_bundled
    return load_bundled(locus, source=source, organism=organism)


def check_junctions(seqs):
    """Split ``seqs`` into (junctions, suspect) by the conserved anchors.

    **CDR3 is not junction.** VDJdb's column is *named* ``cdr3`` but holds junctions — Cys104 and
    Phe118/Trp included. An anchor-stripped IMGT CDR3 scores **exactly 0.0** with no error, so a
    silently mis-typed input reports a precursor frequency of zero rather than failing. Callers
    should check before scoring and report the dropped count rather than letting it vanish.
    """
    ok, bad = [], []
    for s in seqs:
        t = (s or "").strip().upper()
        (ok if t.startswith(_JUNCTION_START) and t.endswith(_JUNCTION_END) else bad).append(t)
    return ok, bad


def pgen(model, junctions, v=None, j=None, mismatches: int = 0, threads: int = 0) -> list[float]:
    """Per-junction `Pgen` — the vector behind every mass in this module.

    ``mismatches=1`` returns the **closed Hamming-1 ball** mass of each junction in closed form
    (inclusion–exclusion inside vdjtools, no enumeration), which is the frequency proxy Pogorelyy
    et al. used. ``v``/``j`` are per-junction allele-resolution call lists or ``None`` to
    marginalise; marginal and conditioned are different quantities, so never mix them in one
    comparison.
    """
    seqs = list(junctions)
    if not seqs:
        return []
    return [float(x) for x in _native().pgen_aa_batch(
        model, seqs, v=v, j=j, mismatches=mismatches, threads=threads)]


def observed_mass(model, junctions, v=None, j=None, threads: int = 0) -> float:
    """Sum of `Pgen` over the given junctions -- a strict lower bound on `F(e)`.

    ``v``/``j`` are per-junction allele-resolution call lists (``TRBV27*01``, not ``TRBV27``), or
    ``None`` to marginalise over V/J. The two are **different quantities** -- the marginal is larger
    -- so never mix them within one comparison.
    """
    return float(sum(pgen(model, junctions, v=v, j=j, threads=threads)))


# --------------------------------------------------------------------------- coverage correction
def _fit_capture_rate(u, m, n):
    """MLE of ``theta`` in ``p_i = 1 - exp(-theta * u_i)`` under a zero-truncated Binomial(n, p_i).

    ``u`` are Pgen values rescaled to unit median so the search bracket is scale-free. Returns
    ``(theta, hit_boundary)``. Grid then golden-section: the likelihood is unimodal in practice but
    not provably so, and a coarse grid makes the refinement immune to a bad initial bracket.
    """
    import numpy as np

    u = np.asarray(u, dtype=float)
    m = np.asarray(m, dtype=float)

    def loglik(x):
        t = np.exp(x) * u
        # log(1 - exp(-t)) via -expm1 keeps precision for small t
        return float(np.sum(m * np.log(-np.expm1(-t)) - (n - m) * t - np.log(-np.expm1(-n * t))))

    grid = np.linspace(math.log(1e-8), math.log(1e8), 321)
    vals = np.array([loglik(x) for x in grid])
    k = int(np.argmax(vals))
    if k == 0 or k == len(grid) - 1:
        return float(np.exp(grid[k])), True

    lo, hi = grid[k - 1], grid[k + 1]
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = hi - phi * (hi - lo), lo + phi * (hi - lo)
    fa, fb = loglik(a), loglik(b)
    for _ in range(60):
        if fa < fb:
            lo, a, fa = a, b, fb
            b = lo + phi * (hi - lo)
            fb = loglik(b)
        else:
            hi, b, fb = b, a, fa
            a = hi - phi * (hi - lo)
            fa = loglik(a)
    return float(math.exp((lo + hi) / 2.0)), False


def coverage_corrected_mass(model, junctions, multiplicity, n_units: int | None = None,
                            threads: int = 0) -> dict:
    """`F(e)` with the size-biased sampling deficit put back — the estimator, not the bound.

    ``multiplicity[i]`` is the number of **independent capture units** (donors, or failing that
    studies) that reported ``junctions[i]`` for this epitope. Take it from donor/study counts, never
    from VDJdb record counts: rows duplicate the same donor's TCR across curation passes, so a
    record count is not a recapture and inflating it collapses the correction. ``n_units`` is the
    total number of capture units for the epitope (defaults to ``max(multiplicity)``).

    **Model.** Each cognate junction is captured by one unit with probability
    ``p_i = 1 - exp(-theta * pi_i)`` — Poisson sampling at a rate proportional to `Pgen`, which is
    exactly the size-biasing that makes ``observed_mass`` biased. ``theta`` is fitted by maximising
    the zero-truncated Binomial(``n_units``, ``p_i``) likelihood over the *observed* junctions
    (unobserved ones contribute nothing, hence the truncation). The corrected mass is then
    Horvitz–Thompson:

    ``F_hat = sum_i pi_i / (1 - (1 - p_i)^n)``

    which is unbiased for the total over the **whole** cognate set, unobserved members included,
    because each observed member is upweighted by its own inclusion probability. Returns
    ``{"observed", "corrected", "coverage", "theta", "gt_coverage", "f1", "n_units", "n_seqs",
    "n_zero", "degenerate", "reason"}``; ``coverage = observed / corrected``.

    **What it assumes.** (1) Capture is independent across units — false where two studies share
    donors or one study's TCRs were re-curated into another, which biases ``theta`` upward and the
    correction downward. (2) The capture curve is the one-parameter saturating form above; real
    capture also depends on assay sensitivity and HLA typing of the cohort, which this folds into a
    single ``theta``. (3) `Pgen` is the right size variable — it is the one the size-bias argument
    names, but clonal expansion in the source samples adds a second, unmodelled one. (4) The cognate
    set is fixed; a junction that is cognate in one donor and not another breaks the frame entirely.

    **Where it breaks, loudly.** Textbook Good–Turing is *known-bad* on TCR data — Laydon et al.,
    *PLoS Comput Biol* 2014;10:e1003646 (PMID 24945836) measure 61.7% median error on real TCR
    abundance data because the capture-probability distribution is far too heterogeneous for the
    uniform-multinomial assumption. That is why capture here is a fitted function of `Pgen` rather
    than the flat ``1 - f1/N``; the flat number is still returned as ``gt_coverage`` so the two can
    be compared, but it is a diagnostic, not the estimate. When **every junction is a singleton**
    there is no recapture information at all, the likelihood is maximised as ``theta -> 0`` and the
    Horvitz–Thompson sum diverges: the function then sets ``degenerate=True``, names the reason and
    returns the ``observed`` bound unchanged rather than an infinity or a ZeroDivisionError. The
    same happens with fewer than two capture units, and if the fit hits the search boundary.
    """
    seqs = list(junctions)
    mult = [int(x) for x in multiplicity]
    if len(mult) != len(seqs):
        raise ValueError(f"multiplicity has {len(mult)} entries for {len(seqs)} junctions")
    if any(x < 1 for x in mult):
        raise ValueError("multiplicity entries must be >= 1: an unobserved junction is not a member")
    n = int(n_units) if n_units is not None else (max(mult) if mult else 0)
    if any(x > n for x in mult):
        raise ValueError(f"multiplicity exceeds n_units={n}: a junction cannot be captured twice "
                         "by the same unit")

    pis = pgen(model, seqs, threads=threads)
    observed = float(sum(pis))
    f1 = sum(1 for x in mult if x == 1)
    captures = sum(mult)
    out = {"observed": observed, "corrected": observed, "coverage": 1.0, "theta": None,
           "gt_coverage": (1.0 - f1 / captures) if captures else 0.0,
           "f1": f1, "n_units": n, "n_seqs": len(seqs), "n_zero": sum(1 for p in pis if p <= 0),
           "degenerate": True, "reason": None}

    if not seqs:
        out["reason"] = "no junctions"
        return out
    if n < 2:
        out["reason"] = "fewer than two capture units: no recapture information"
        return out
    if f1 == len(mult):
        out["reason"] = ("every junction is a singleton: the capture curve is unidentified "
                         "(theta -> 0, the Horvitz-Thompson sum diverges)")
        return out

    keep = [i for i, p in enumerate(pis) if p > 0]
    if not keep:
        out["reason"] = "every junction has Pgen 0 under this model"
        return out
    kept = [pis[i] for i in keep]
    med = sorted(kept)[len(kept) // 2]
    u = [p / med for p in kept]
    theta, boundary = _fit_capture_rate(u, [mult[i] for i in keep], n)
    if boundary:
        out["reason"] = "the capture-rate fit hit the search boundary; theta is not identified"
        return out

    corrected = 0.0
    for p, ui in zip(kept, u):
        incl = -math.expm1(-n * theta * ui)
        if incl <= 0:                                # pragma: no cover - guarded by `boundary`
            out["reason"] = "an inclusion probability underflowed to 0"
            return out
        corrected += p / incl
    out.update(corrected=corrected, coverage=observed / corrected if corrected > 0 else 1.0,
               theta=theta / med, degenerate=False)
    return out


# --------------------------------------------------------------------------- neighbourhood mass
def ball_mass(model, junctions, r: int = 1, threads: int = 0) -> dict:
    """Mass of the **union** of Hamming-``r`` balls around ``junctions``.

    Returns ``{"union", "naive_sum", "overlap", "n_union", "n_seqs"}`` where ``naive_sum`` adds the
    per-sequence closed-ball masses independently and ``overlap = 1 - union/naive_sum`` is the
    fraction that double-counting would have invented.

    V/J are deliberately **not** accepted: a substituted neighbour need not keep the centre's V/J
    assignment, so conditioning the ball on the centre's call would be wrong. This marginalises.
    """
    from seqtree.distance import neighbourhood_union
    native = _native()
    seqs = [s for s in junctions if s]
    if not seqs:
        return {"union": 0.0, "naive_sum": 0.0, "overlap": 0.0, "n_union": 0, "n_seqs": 0}

    members = neighbourhood_union(seqs, r=r)
    union = float(sum(native.pgen_aa_batch(model, list(members), threads=threads)))

    if r == 1:                                   # closed-ball closed form, no enumeration
        naive = float(sum(native.pgen_aa_batch(model, seqs, mismatches=1, threads=threads)))
    else:
        naive = 0.0
        for s in seqs:
            naive += float(sum(native.pgen_aa_batch(
                model, list(neighbourhood_union([s], r=r)), threads=threads)))
    return {"union": union, "naive_sum": naive,
            "overlap": (1.0 - union / naive) if naive > 0 else 0.0,
            "n_union": len(members), "n_seqs": len(seqs)}


def shell_profile(model, junctions, r: int = 1, alpha: float = ALPHA_PER_EDIT,
                  threads: int = 0, max_members: int = MAX_BALL_MEMBERS) -> dict:
    """`ball_mass` resolved by exact edit distance, with cognacy retention applied per shell.

    A ball at radius ``r`` treats a junction ``r`` substitutions from an observed cognate TCR as
    fully cognate, which it is not. Shell ``k`` is the set of sequences whose distance to the
    **nearest** observed junction is exactly ``k`` (seqtree's ``shell=`` semantics — min-distance,
    so the shells partition the ball), and the retained estimate is

    ``F ≈ sum_k alpha**k * mass(shell k)``

    with ``alpha`` the per-edit cognacy retention, default :data:`ALPHA_PER_EDIT` (0.1, Mayer &
    Callan 2023). ``alpha=1`` reproduces the raw union; ``alpha=0`` collapses to
    :func:`observed_mass`.

    Returns ``{"shells": [{"r", "n", "mass", "alpha"}...], "retained", "union", "n_union",
    "n_seqs", "alpha"}``.

    **Memory.** The shells are materialised as Python strings. ``r=1`` costs ~19L sequences per
    junction and is always fine; ``r=2`` is ~180*L^2/2 each and is not — 300 measured junctions at
    ``r=2`` come to ~9.9M sequences / ~1.8 GB. The union is sized with seqtree's ``union_size``
    *before* anything is enumerated, and a request above ``max_members``
    (:data:`MAX_BALL_MEMBERS`) raises ``MemoryError``. Split the junction list and add the shell
    masses if you need a bigger radius; shells are disjoint, so that is exact.
    """
    from seqtree.distance import neighbourhood_union, union_size
    native = _native()
    seqs = [s for s in junctions if s]
    if not seqs:
        return {"shells": [], "retained": 0.0, "union": 0.0, "n_union": 0, "n_seqs": 0,
                "alpha": alpha}
    if r < 0:
        raise ValueError("r must be >= 0")

    n_union = union_size(seqs, r=r)
    if n_union > max_members:
        raise MemoryError(
            f"the radius-{r} union of {len(seqs)} junctions holds {n_union:,} sequences, above "
            f"max_members={max_members:,} (~{n_union * 190 / 1e9:.1f} GB as Python strings). "
            "Lower r, split the junction list and add the shell masses, or raise max_members.")

    shells, retained = [], 0.0
    for k in range(r + 1):
        members = neighbourhood_union(seqs, r=k, shell=True) if k else list(dict.fromkeys(seqs))
        mass = float(sum(native.pgen_aa_batch(model, list(members), threads=threads)))
        w = alpha ** k
        shells.append({"r": k, "n": len(members), "mass": mass, "alpha": w})
        retained += w * mass
    return {"shells": shells, "retained": retained,
            "union": float(sum(s["mass"] for s in shells)),
            "n_union": n_union, "n_seqs": len(seqs), "alpha": alpha}


# --------------------------------------------------------------------------- cluster PWM motifs
@dataclass(frozen=True)
class ClusterMotif:
    """One VDJdb cluster PWM as a degenerate motif ready for :func:`motif_mass`.

    ``allowed`` is one string of permitted residues per position (``""`` = wildcard). The cluster is
    V/J/length-pinned, so ``v``/``j`` are the conditioning to pass alongside it.
    """
    cid: str
    epitope: str
    gene: str
    species: str
    v: str
    j: str
    length: int
    size: int
    allowed: tuple[str, ...]


def load_cluster_motifs(path, threshold: float = MOTIF_FREQ_THRESHOLD, species: str | None = None,
                        gene: str | None = None, epitope: str | None = None,
                        min_size: int = 2) -> list[ClusterMotif]:
    """Read VDJdb's ``motif_pwms.txt`` into per-position allowed-residue sets.

    The file is one row per (cluster, position, residue) with ``freq`` the residue's frequency
    *within* the cluster; only residues actually observed are listed. Because a cluster is pinned to
    one V, one J and one length, thresholding ``freq`` per position yields exactly the ``allowed``
    argument :func:`motif_mass` wants — no alignment, no register search, no enumeration.

    ``threshold`` is the per-position frequency cut, default :data:`MOTIF_FREQ_THRESHOLD`. A
    position never comes back empty: if no residue clears the cut the modal residue(s) are kept, so
    raising the threshold shrinks the motif monotonically towards the consensus sequence rather than
    zeroing its mass. Positions absent from the file become wildcards.

    ``species``/``gene``/``epitope`` filter exactly; ``min_size`` drops clusters below that many
    members (``csz``). Accepts a plain or gzipped path.
    """
    p = Path(path)
    opener = (lambda: gzip.open(p, "rt")) if p.suffix == ".gz" else (lambda: open(p, "r"))
    rows: dict[str, dict] = {}
    with opener() as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            cid = rec["cid"]
            c = rows.get(cid)
            if c is None:
                c = rows[cid] = {"meta": rec, "pos": {}}
            c["pos"].setdefault(int(rec["pos"]), {})[rec["aa"]] = float(rec["freq"])

    out = []
    for cid, c in rows.items():
        meta = c["meta"]
        if species is not None and meta["species"] != species:
            continue
        if gene is not None and meta["gene"] != gene:
            continue
        if epitope is not None and meta["antigen.epitope"] != epitope:
            continue
        size = int(meta["csz"])
        if size < min_size:
            continue
        length = int(meta["len"])
        allowed = []
        for i in range(length):
            d = c["pos"].get(i)
            if not d:
                allowed.append("")                          # unobserved position -> wildcard
                continue
            keep = sorted(a for a, f in d.items() if f >= threshold)
            if not keep:                                    # never emit an empty set
                top = max(d.values())
                keep = sorted(a for a, f in d.items() if f >= top)
            allowed.append("".join(keep))
        out.append(ClusterMotif(
            cid=cid, epitope=meta["antigen.epitope"], gene=meta["gene"], species=meta["species"],
            v=meta["v.segm.repr"], j=meta["j.segm.repr"], length=length, size=size,
            allowed=tuple(allowed)))
    return out


def motif_mass(model, allowed, v=None, j=None) -> float:
    """`Pgen` of every junction matching a degenerate motif.

    ``allowed`` is one entry per position, each a string of permitted residues; ``""`` or ``"X"``
    means any residue. A VDJdb cluster PWM is V/J/length-pinned, so thresholding it per position
    gives exactly this — and one call returns the whole cluster's mass with no enumeration and no
    inclusion–exclusion.

    Unlike :func:`observed_mass` and :func:`ball_mass` this scores a **set**, so it carries no
    observed-sample coverage bias. Pass the motif's own ``v``/``j`` — cluster motifs are V/J-pinned
    and the conditioned quantity is the right one for them.
    """
    return float(_native().pgen_aa_degenerate(model, list(allowed), v=v, j=j))


# --------------------------------------------------------------------------- the A-vs-B check
def cross_check(model, junctions, allowed, r: int = 1, alpha: float = ALPHA_PER_EDIT,
                threads: int = 0, max_members: int = MAX_BALL_MEMBERS) -> dict:
    """Two independent estimates of the same `F(e)` — and their disagreement is the missing mass.

    Route **A** (:func:`motif_mass`) scores a *set*: it asks the recombination model for the total
    mass of every junction the motif admits, so it never touches the observed sample and carries
    none of its coverage bias. Route **B** (:func:`observed_mass` / :func:`shell_profile`) starts
    from the junctions VDJdb actually recorded, so it is bounded below by how deeply the epitope was
    sampled. Both are estimates of `F(e)` for the same epitope.

    Returns ``{"set_mass", "observed_mass", "ball_mass", "retained_mass", "ratio_observed",
    "ratio_retained", "missing_fraction", "n_seqs", "n_union"}``.

    **Interpretation.** ``ratio_observed = set_mass / observed_mass`` is the factor by which the
    sample under-counts, and ``missing_fraction = 1 - observed_mass/set_mass`` the share of the
    cognate mass never observed. ``ratio_retained`` repeats it against the shell-weighted estimate:
    if the neighbourhood correction is doing its job, ``ratio_retained`` is closer to 1 than
    ``ratio_observed`` — that is the whole claim of the ``ball``/``shell`` route, tested rather than
    assumed. A ratio **below 1** is informative in the other direction: the motif is tighter than
    the sample it was built from, i.e. the threshold is too strict or the cluster is a proper subset
    of the epitope's cognate TCRs (which it usually is — one epitope has several clusters).

    Both routes are **marginalised over V/J** so the two numbers are the same quantity. Cluster
    motifs are V/J-pinned and :func:`motif_mass` will condition on request, but a conditioned A
    against a marginal B is a category error and is not offered here.
    """
    set_mass = motif_mass(model, allowed)
    obs = observed_mass(model, junctions, threads=threads)
    prof = shell_profile(model, junctions, r=r, alpha=alpha, threads=threads,
                         max_members=max_members)
    return {"set_mass": set_mass, "observed_mass": obs, "ball_mass": prof["union"],
            "retained_mass": prof["retained"],
            "ratio_observed": (set_mass / obs) if obs > 0 else float("inf"),
            "ratio_retained": (set_mass / prof["retained"]) if prof["retained"] > 0 else float("inf"),
            "missing_fraction": (1.0 - obs / set_mass) if set_mass > 0 else 0.0,
            "n_seqs": prof["n_seqs"], "n_union": prof["n_union"]}


def demo() -> None:
    """Self-check: ``python -m mhcmatch.precursor`` (skips cleanly without vdjtools)."""
    try:
        _native()
    except ImportError as e:
        print(f"skip - {e}")
        return

    m = load_model("TRB")
    a = "CASSLAPGATNEKLFF"
    b = "CASSLAPGATNEKLYF"                       # Hamming-1 from a
    far = "CASSQDRDTQYF"                         # different length -- balls cannot overlap

    ok, bad = check_junctions([a, "ASSLAPGATNEKL", b])
    assert ok == [a, b] and bad == ["ASSLAPGATNEKL"], (ok, bad)

    # the lower bound is a plain sum
    s = observed_mass(m, [a, b])
    assert s > 0 and abs(s - (observed_mass(m, [a]) + observed_mass(m, [b]))) < 1e-30

    # a ball is bigger than its centre, and the union is smaller than the naive sum when the
    # centres are close enough for their balls to intersect
    near = ball_mass(m, [a, b])
    assert near["union"] > s, "the ball must exceed the observed mass"
    assert near["union"] < near["naive_sum"], near
    assert near["overlap"] > 0, near

    # ... and there is nothing to double-count when the balls are disjoint by construction
    disj = ball_mass(m, [a, far])
    assert abs(disj["overlap"]) < 1e-12, disj

    # a single sequence: the union IS the closed ball, so the two agree exactly
    one = ball_mass(m, [a])
    assert abs(one["union"] - one["naive_sum"]) / one["naive_sum"] < 1e-9, one

    # the shells partition the ball, and retention interpolates between the bound and the union
    prof = shell_profile(m, [a, b], r=1)
    assert abs(prof["union"] - near["union"]) / near["union"] < 1e-9, (prof, near)
    assert s < prof["retained"] < prof["union"], prof

    # motif mass: pinning every position reproduces the single-sequence Pgen
    assert abs(motif_mass(m, list(a)) - observed_mass(m, [a])) < 1e-30
    # widening one position can only add mass
    wide = list(a)
    wide[5] = ""
    assert motif_mass(m, wide) > motif_mass(m, list(a))

    # ... and the A-vs-B ratio is exactly 1 when the "sample" is the whole set
    import itertools
    allowed = ["C", "A", "S", "S", "LM", "AG", "N", "E", "K", "L", "F", "F"]
    members = ["".join(c) for c in itertools.product(*allowed)]
    xc = cross_check(m, members, allowed)
    assert abs(xc["ratio_observed"] - 1.0) < 1e-9, xc

    print(f"ok - observed {s:.3e} | union {near['union']:.3e} over {near['n_union']} seqs "
          f"| double-counting avoided {near['overlap']:.1%} "
          f"| retained (alpha={ALPHA_PER_EDIT}) {prof['retained']:.3e}")


if __name__ == "__main__":
    demo()
