"""T-cell precursor frequency for an epitope: how much repertoire mass can see it.

The estimand is

    F(e) = sum over the cognate set C_e of pi(tau)

— the probability that a random naive-repertoire junction recognises epitope ``e``. This is the
continuous quantity behind "immunogenic": recognition is a spectrum, and `F(e)` is where on it a
given epitope sits.

Three estimators of the same `F(e)`, in increasing order of model commitment:

``observed_mass``
    Sum of `Pgen` over the junctions actually recorded for the epitope. A **strict lower bound**,
    and a biased one: VDJdb samples cognate TCRs **size-biased by Pgen** (a TCR enters the record
    roughly in proportion to its repertoire frequency), so the observed members are systematically
    the high-`Pgen` ones and the deficit does not shrink with more studies at the same depth.

``ball_mass``
    Mass of the **union** of Hamming-`r` balls around the observed junctions. Cognate TCRs are
    near-duplicates by construction, so their balls overlap and the *sum* of per-sequence ball
    masses double-counts; the union is the correct object. The returned ``overlap`` quantifies that
    double-counting directly, which is worth reporting rather than hiding — it is a measurement of
    how tight the specificity group is.

``motif_mass``
    `Pgen` of a degenerate motif (a VDJdb cluster PWM is V/J/length-pinned, hence exactly a
    per-position residue set). This computes the mass of a **set** directly, so unlike the two above
    it does not suffer the observed-sample coverage bias at all.

That last point is the useful one: ``motif_mass`` and ``ball_mass`` estimate the same quantity by
independent routes, so **their disagreement is an estimate of the missing mass**.

Requires the optional ``vdjtools`` dependency (the recombination model). Nothing here reimplements
`Pgen` — the DP, the closed Hamming-1 ball and the degenerate/masked DP all live in vdjtools, and
the neighbourhood enumeration lives in seqtree.
"""

from __future__ import annotations

__all__ = ["load_model", "check_junctions", "observed_mass", "ball_mass", "motif_mass"]

_JUNCTION_START = "C"
_JUNCTION_END = ("F", "W")


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


def observed_mass(model, junctions, v=None, j=None, threads: int = 0) -> float:
    """Sum of `Pgen` over the given junctions -- a strict lower bound on `F(e)`.

    ``v``/``j`` are per-junction allele-resolution call lists (``TRBV27*01``, not ``TRBV27``), or
    ``None`` to marginalise over V/J. The two are **different quantities** -- the marginal is larger
    -- so never mix them within one comparison.
    """
    native = _native()
    seqs = list(junctions)
    if not seqs:
        return 0.0
    return float(sum(native.pgen_aa_batch(model, seqs, v=v, j=j, threads=threads)))


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

    # motif mass: pinning every position reproduces the single-sequence Pgen
    assert abs(motif_mass(m, list(a)) - observed_mass(m, [a])) < 1e-30
    # widening one position can only add mass
    wide = list(a)
    wide[5] = ""
    assert motif_mass(m, wide) > motif_mass(m, list(a))

    print(f"ok - observed {s:.3e} | union {near['union']:.3e} over {near['n_union']} seqs "
          f"| double-counting avoided {near['overlap']:.1%}")


if __name__ == "__main__":
    demo()
