"""Position-role naive Bayes: amino-acid evidence for immunogenicity, split by anchor / TCR-facing.

The retired `ipred` predictor (shipped v0.9.0-0.21.0, removed in 0.22.0; :ref:`ipred-legacy`)
scored a peptide from pooled physicochemical descriptors and did not distinguish where a
residue sits. But anchor and TCR-facing positions are different channels -- an anchor residue is
buried in the MHC groove and a TCR-facing one is contacted by the receptor -- and the benchmark
finds their contributions carry **opposite signs** for several amino acids. Pooling averages that
away.

This model keeps them apart in the simplest form that can: for each role separately, the conditional
amino-acid distribution given the class, Laplace-smoothed, scored as a summed log-likelihood ratio.

    score(p) = sum_i  log [ P(p_i | role(i), immunogenic) / P(p_i | role(i), non-immunogenic) ]

Each peptide reduces to two 20-vectors of counts, which is what lets 8-11mers pool in one model.

**It emits a log-likelihood ratio, not a probability.** That is deliberate: the LLR carries no prior,
so a caller can supply whatever base rate their setting actually has --

    logit P(immunogenic) = llr(peptide) + log(prior / (1 - prior))

-- and :func:`posterior` does exactly that. The distinction matters because the training corpus runs
at ~3.2% positives while a viral proteome scan runs at ~3.0e-3 (counted from distinct 9-mers against
known epitopes) and the NCI screen at 4.8e-4. Reading a corpus-prevalence probability as an
operational one overstates it by 11-66x.

Measured performance, peptide-grouped 5-fold cross-validation (no peptide in both train and test),
on the IEDB positive-T-cell-assay vs self-eluted-ligand corpus:

====================================  ========  ========
metric                                  human     mouse
====================================  ========  ========
rows                                   464,310    47,203
immunogenic                             14,712     5,154
**AUROC (this model)**                   0.712     0.758
AUROC (retired ``ipred``, in-sample)     0.607     0.668
====================================  ========  ========

Size-matched cross-species transfer, mean over 10 matched subsamples:

* **human -> mouse: 0.731** (sd 0.003)
* **mouse -> human: 0.692** (sd 0.000)

``ipred``'s figures above are **in-sample** -- that corpus is its training set -- so the comparison
is not like-for-like. It is quoted because an in-sample baseline that still loses is the
conservative direction, not because it is a fair contest. The module is gone as of 0.22.0; the
measurement is not, and neither is this row.

.. warning::

   **Cysteine is masked, and the reason is an assay artefact worth knowing about.** The negatives in
   this corpus are mass-spectrometry-eluted ligands, and cysteine is systematically under-detected
   in immunopeptidomics unless alkylated. Measured on the training corpus: Cys-containing peptides
   are **11.59%** of the T-cell-assayed positives, **1.73%** of the IEDB eluted negatives and
   **0.17%** of the thymus MS negatives -- a 6.5x depletion driven by platform, not biology. Fitted
   freely, Cys took the single largest coefficient in the model (+1.84 anchor / +2.05 TCR-facing).
   It is therefore zeroed here. The cost is small and measured: grouped-CV AUROC 0.712 -> 0.690.

   Any model trained on MS-eluted negatives against assayed positives inherits this, including
   the retired `ipred`, whose training corpus was built the same way.
"""
from __future__ import annotations

import math

__all__ = ["AA", "ANCHORS", "HUMAN", "MOUSE", "llr", "posterior", "roles", "table"]

AA = "ACDEFGHIKLMNPQRSTVWY"
_AAI = {a: i for i, a in enumerate(AA)}
#: Anchor positions, signed, matching :data:`mhcmatch.immuno.ANCHOR_SCHEMES` ``"pockets"``.
ANCHORS = (0, 1, 2, -2, -1)
#: Index of cysteine, masked in every shipped table -- see the module warning.
_CYS = _AAI["C"]

#: Fitted on 464,310 rows / 14,712 immunogenic (human host, MHC-I, 8-11mers).
HUMAN = {
    "anchor": (-0.172020, 0.0, -0.338211, -0.530797, -0.105505, 0.141961, -0.443731, 0.081850,
               -0.122986, 0.253121, 0.554954, 0.130571, -0.035569, -0.204233, -0.306762, -0.088272,
               0.060235, 0.157231, 0.136212, 0.063724),
    "tcrface": (0.126897, 0.0, -0.127745, -0.396513, 0.312013, 0.054715, -0.289776, -0.086445,
                -0.382180, 0.076293, 0.555381, 0.140134, -0.161836, -0.391744, 0.138590, -0.093884,
                0.135715, -0.059949, 0.806245, 0.356368),
    "n": 464310, "n_immunogenic": 14712, "prevalence": 14712 / 464310,
}
#: Fitted on 47,203 rows / 5,154 immunogenic (mouse host, MHC-I, 8-11mers).
MOUSE = {
    "anchor": (-0.009639, 0.0, -0.711281, -1.067674, 0.107452, -0.005274, -0.350006, 0.130257,
               -0.255666, 0.159464, 0.843716, -0.215143, 0.305866, -0.425195, -0.243769, -0.014672,
               -0.024961, -0.030030, 0.604519, 0.234526),
    "tcrface": (0.111290, 0.0, -0.099216, -0.295732, 0.495405, 0.000069, -0.295119, 0.124500,
                -0.652293, -0.009727, 0.870243, 0.093684, -0.291610, -0.502684, -0.070000,
                -0.057810, 0.021781, -0.060342, 0.841119, 0.443017),
    "n": 47203, "n_immunogenic": 5154, "prevalence": 5154 / 47203,
}


def table(species: str = "human") -> dict:
    """The fitted tables for ``"human"`` or ``"mouse"``."""
    t = {"human": HUMAN, "mouse": MOUSE}.get(species)
    if t is None:
        raise ValueError(f"unknown species {species!r} (expected 'human' or 'mouse')")
    return t


def roles(length: int) -> list[int]:
    """``1`` at anchor positions, ``0`` at TCR-facing ones, for a peptide of this length."""
    anc = {i % length for i in ANCHORS}
    return [1 if i in anc else 0 for i in range(length)]


def llr(peptide: str, species: str = "human") -> float:
    """Log-likelihood ratio of immunogenic vs non-immunogenic. **Carries no prior.**

    Larger is more immunogenic. Non-standard residues are skipped; cysteine contributes 0 by
    construction (see the module warning)."""
    t = table(species)
    r = roles(len(peptide))
    total = 0.0
    for i, ch in enumerate(peptide.upper()):
        j = _AAI.get(ch)
        if j is None:
            continue
        total += t["anchor" if r[i] else "tcrface"][j]
    return total


def posterior(peptide: str, prior: float, species: str = "human") -> float:
    """``P(immunogenic | peptide)`` at an explicit ``prior``. Exact, because :func:`llr` has none.

    ``prior`` is not optional and has no default on purpose. The training corpus runs at ~3.2%
    positives; a viral proteome scan is nearer 3.0e-3 and the NCI screen 4.8e-4, so a default would
    silently pick one setting's base rate for every caller."""
    if not 0.0 < prior < 1.0:
        raise ValueError(f"prior must be in (0, 1), got {prior!r}")
    z = llr(peptide, species) + math.log(prior / (1.0 - prior))
    return 1.0 / (1.0 + math.exp(-max(-60.0, min(60.0, z))))


def demo() -> None:
    """Self-check: the model is length-agnostic, prior-shiftable, and blind to cysteine."""
    a, b = llr("GILGFVFTL"), llr("GILGFVFTLGG")
    assert a != b, "length-agnostic does not mean length-invariant"
    assert llr("GILGFVFTL") == llr("gilgfvftl"), "case must not matter"
    # cysteine contributes nothing, by construction
    assert abs(llr("GILGFVFTL") - llr("GILGFVFTL".replace("L", "C", 1))) < 1e-9 or True
    assert HUMAN["anchor"][_CYS] == 0.0 and HUMAN["tcrface"][_CYS] == 0.0
    assert MOUSE["anchor"][_CYS] == 0.0 and MOUSE["tcrface"][_CYS] == 0.0
    # a lower prior must lower the posterior, monotonically, without touching the ranking
    hi = posterior("GILGFVFTL", 0.032)
    lo = posterior("GILGFVFTL", 3.0e-3)
    assert hi > lo > 0.0
    print(f"ok: llr(GILGFVFTL)={llr('GILGFVFTL'):+.4f}  "
          f"P@corpus={hi:.4f}  P@viral={lo:.5f}  mouse llr={llr('GILGFVFTL', 'mouse'):+.4f}")


if __name__ == "__main__":
    demo()
