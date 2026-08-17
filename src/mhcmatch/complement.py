"""Complementarity: how well a presented peptide complements a T-cell repertoire.

This is the recognition axis, and it is what :mod:`mhcmatch.ipred` was reaching for. `ipred` summed
two principal components of the amino-acid property matrix over the whole peptide and added its
length -- 13 parameters, fitted by EM. That construction is kept and five blocks are added, each
answering something the pooled version provably cannot:

============  =========================================================================
block         what it adds
============  =========================================================================
``phys``      PC1/PC2 of the property matrix summed over the peptide, plus length. **The
              `ipred` feature set**, kept as the floor everything else is measured against.
``role``      The same components computed **separately over MHC-facing and TCR-facing**
              residues, plus Kidera KF4 (hydropathy) per role. The two channels carry
              opposite-sign contributions for several residues; a pooled sum reports
              their difference, weighted by the corpus's composition.
``pot``       Contact potentials, one per side: **MJ1996** on the anchors (burial in a
              groove -- MJ is 96.4% one-body with hydropathy as its dominant mode) and
              **TCRen marginalised over a real CDR3 repertoire** on the TCR-facing
              residues. TCRen is only 3.29% one-body, so no per-residue scale can be
              extracted from it and the unknown receptor side is integrated out instead.
``motif``     Contiguity: longest run, number of runs and above-median fraction of
              hydrophobic TCR-facing residues. A run of 3-4 is a different object from
              the same residues scattered, and **no sum can express the difference**.
``aa``        Residue **identity** as a log-odds per amino acid per role. Every block
              above projects onto a property; this one does not. Its two columns sum to
              exactly :func:`mhcmatch.posbayes.llr`.
``kmer``      The same over **adjacent TCR-facing residue pairs** -- a preference for a
              specific dipeptide that no marginal composition feature can express.
============  =========================================================================

**The head is linear, and that is a measured choice.** The shipped ``posbayes`` score is a *sum* of
the two role log-odds -- weights fixed at 1 on two of these columns. A diagonal-covariance Gaussian
classifier cannot represent that: it maps each column through its own quadratic and re-weights by
inverse class variances, so the additive form is outside its hypothesis space and the extra blocks
are paid for out of a worse fit to the term carrying most of the signal. On the training corpus the
EM Gaussian reaches 0.657 grouped-CV AUROC on the ``aa`` block where the plain sum reaches 0.711.
A linear head *contains* the sum as a special case, so whatever it adds is genuinely an addition.
The EM Gaussian parameters are vendored alongside anyway, so the comparison stays re-checkable.

**It emits a log-odds, not a probability.** Like :mod:`~mhcmatch.posbayes`, the corpus's own base
rate is divided out, so a caller supplies whatever prevalence their setting actually has::

    logit P(immunogenic) = score(peptide) + log(prior / (1 - prior))

and :func:`posterior` does exactly that.

**Vectorised.** :func:`score` takes an iterable and returns an array; the whole feature set is two
(n, 20) count matrices times a handful of property vectors, so scoring a full deposit is seconds,
not minutes. Pass a list, not a loop.

    >>> from mhcmatch import complement
    >>> complement.score(["GILGFVFTL", "SIINFEKL"])            # doctest: +SKIP
    array([1.79, 0.42])

Parameters are vendored in ``mhcmatch/data/complement_mhc1.json`` and never refitted at import.
Provenance, the block-by-block cross-validation, the corpus-transfer matrix and the size-matched
cross-species transfer are in the benchmark repo (``bench/results/complementarity.md``).

.. warning::

   **Class I only.** The role split is the class-I one (P1-P3, PΩ-1, PΩ). A class-II ligand is
   anchored by the P1/P4/P6/P9 core of a 9-mer register floating inside a longer peptide
   (:func:`mhcmatch.store.anchor_indices`), so applying this scheme to it labels the wrong residues
   as anchors and returns a confident, wrong number.
"""
from __future__ import annotations

import json
import math
from importlib import resources
from statistics import median

import numpy as np

from . import ipred
from .data import aa_tables

__all__ = ["AA", "ANCHORS", "PARATOPE", "PARAMS", "BLOCKS", "encode", "feature_names",
           "features", "score", "posterior", "parameters"]

AA = "ACDEFGHIKLMNPQRSTVWY"
_AAI = {a: i for i, a in enumerate(AA)}
#: MHC-facing positions, signed, matching :data:`mhcmatch.immuno.ANCHOR_SCHEMES` ``"pockets"`` and
#: :data:`mhcmatch.posbayes.ANCHORS`.
ANCHORS = (0, 1, 2, -2, -1)

#: TCRen marginalised over 28,250,990 TRB IMGT CDR3 loops -- ``paratope(a) = sum_b f(b) TCRen(b,a)``
#: -- and its spread over the same distribution. A residue can have a mild mean energy and still be
#: highly discriminating across receptors, so both are features. Measured in the benchmark repo
#: (``bench/results/paratope_basis.md``); the 32M-clonotype repertoire is not needed at runtime.
PARATOPE = {
    "L": (-0.0251, 0.1864), "G": (-0.0000, 0.2369), "S": (0.0124, 0.2607), "V": (0.0143, 0.2991),
    "Y": (0.0172, 0.2474), "M": (0.0345, 0.3921), "Q": (0.0384, 0.3411), "A": (0.0468, 0.3668),
    "P": (0.0469, 0.2501), "I": (0.0486, 0.4160), "R": (0.0515, 0.4210), "F": (0.0754, 0.3616),
    "W": (0.0847, 0.4427), "E": (0.1014, 0.3130), "N": (0.1193, 0.3291), "K": (0.1234, 0.4516),
    "T": (0.1307, 0.3778), "D": (0.1813, 0.5419), "C": (0.1994, 0.5024), "H": (0.2189, 0.5621),
}

#: Feature blocks, in the order the benchmark's cumulative ablation adds them.
BLOCKS = {
    "phys": ["pc1", "pc2", "length"],
    "role": ["pc1_anchor", "pc2_anchor", "pc1_tcr", "pc2_tcr", "kf4_anchor", "kf4_tcr"],
    "pot": ["mj_anchor", "mj_tcr", "para_tcr", "para_sd_tcr"],
    "motif": ["kd_run_max", "kd_run_n", "kd_run_frac"],
    "aa": ["aa_anchor", "aa_tcr"],
    "kmer": ["kmer_llr"],
}
#: Columns computed from a fitted log-odds table rather than from the peptide alone:
#: ``feature -> which count matrix it weights``.
FITTED = {"aa_anchor": "anchor", "aa_tcr": "tcr", "kmer_llr": "pair"}

_SRC = "complement_mhc1.json"


def _load() -> dict:
    with resources.files("mhcmatch.data").joinpath(_SRC).open() as fh:
        p = json.load(fh)
    k = len(p["features"])
    if len(p["logistic"]["coef"]) != k:
        raise ValueError(f"{_SRC}: {len(p['logistic']['coef'])} coefficients for {k} features")
    if len(p["standardizer"]["mean"]) != k or len(p["standardizer"]["std"]) != k:
        raise ValueError(f"{_SRC}: standardizer does not cover {k} features")
    if not 0.0 < p["prevalence"] < 1.0:
        raise ValueError(f"{_SRC}: prevalence {p['prevalence']!r} is not a base rate")
    for name, src in p["log_odds_source"].items():
        want = {"anchor": 20, "tcr": 20, "pair": 400}[src]
        if len(p["log_odds"][name]) != want:
            raise ValueError(f"{_SRC}: {name} has {len(p['log_odds'][name])} cells, want {want}")
    return p


#: The frozen model: standardizer, linear head, the two fitted log-odds tables, and the EM /
#: supervised Gaussian parameters kept for comparison.
PARAMS: dict = _load()


def _scale_vec(tab: dict) -> np.ndarray:
    return np.array([tab[a] for a in AA], dtype=float)


def _basis() -> dict[str, np.ndarray]:
    rs = ipred.residue_scores()
    return {
        "pc1": np.array([rs[a][0] for a in AA]),
        "pc2": np.array([rs[a][1] for a in AA]),
        "kf4": _scale_vec(aa_tables.DESCRIPTORS["KIDERA"]["KF4"]),
        "mj": _scale_vec(aa_tables.MJ_PARTITION),
        "para": np.array([PARATOPE[a][0] for a in AA]),
        "para_sd": np.array([PARATOPE[a][1] for a in AA]),
        "kd": _scale_vec(aa_tables.HYDROPHOBICITY["KyteDoolittle"]),
    }


BASIS = _basis()
#: The "hydrophobic" cut for the run features: the median of the Kyte-Doolittle scale itself, so it
#: is a property of the scale rather than a tuned constant. Same rule as
#: :func:`mhcmatch.immuno._aggregate`.
#:
#: Taken from the plain dict with the stdlib rather than from :data:`BASIS` with ``numpy``, because
#: sphinx mocks ``numpy`` at doc-build time and a module-level ``float(np.median(...))`` then raises
#: on a Mock -- the whole module fails to import and its page renders empty.
KD_THRESHOLD = float(median(aa_tables.HYDROPHOBICITY["KyteDoolittle"].values()))


def feature_names() -> list[str]:
    """Column order of the design matrix, matching the vendored coefficients."""
    return list(PARAMS["features"])


def parameters() -> dict:
    """The fitted model as a plain dict (a copy of the vendored file)."""
    return json.loads(json.dumps(PARAMS))


def encode(peptides) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """``(features, counts)`` for an iterable of peptides -- the vectorised feature builder.

    Everything additive is a matrix product against a residue vector, so the peptide-side feature
    set is two (n, 20) count matrices (one per role) times seven property vectors. Run statistics
    need position order and use a loop over **positions** (at most 11 iterations), never over
    peptides. Peptides are grouped by length so each group is one batch of array operations.

    ``counts`` carries ``anchor`` and ``tcr`` (n, 20) -- the matrices the ``aa`` log-odds tables
    weight -- plus the adjacent TCR-facing residue pairs as a **sparse pair list**, ``pair_code``
    (which of 400 pairs) and ``pair_row`` (which peptide). A 9-mer has at most 3 such pairs, so a
    dense (n, 400) matrix would be 99% zeros and would cost 1.5 GB of temporaries per pass on a
    500k-peptide corpus; :func:`apply_log_odds` sums the sparse form instead.

    Non-standard residues (``X`` masks, ``B``/``J``/``O``/``U``/``Z``) contribute to no count and
    break no run, but still count toward ``length``.
    """
    peps = [str(p).strip().upper() for p in peptides]
    n = len(peps)
    c_anc = np.zeros((n, 20))
    c_tcr = np.zeros((n, 20))
    pair_code: list = []
    pair_row: list = []
    out = {k: np.zeros(n) for k in ("kd_run_max", "kd_run_n", "kd_run_frac")}
    out["length"] = np.array([len(p) for p in peps], dtype=float)

    by_len: dict[int, list[int]] = {}
    for i, p in enumerate(peps):
        by_len.setdefault(len(p), []).append(i)

    for L, idx in by_len.items():
        if L < 3:
            continue
        rows = np.array(idx)
        code = np.array([[_AAI.get(ch, -1) for ch in peps[i]] for i in rows], dtype=np.int16)
        anc = np.zeros(L, dtype=bool)
        for a in ANCHORS:
            anc[a % L] = True
        valid = code >= 0
        safe = np.where(valid, code, 0)

        for mask, mat in ((anc, c_anc), (~anc, c_tcr)):
            sel = valid & mask[None, :]
            if sel.any():
                np.add.at(mat, (np.repeat(rows, L)[sel.ravel()], safe[sel]), 1)

        # Hydrophobic runs, TCR-facing only. A buried anchor BREAKS a run rather than bridging it:
        # from the receptor's point of view two stretches either side of a buried residue are two
        # stretches, not one.
        above = valid & (~anc)[None, :] & (BASIS["kd"][safe] > KD_THRESHOLD)
        cur = np.zeros(len(rows))
        best = np.zeros(len(rows))
        nrun = np.zeros(len(rows))
        prev = np.zeros(len(rows), dtype=bool)
        for j in range(L):
            a = above[:, j]
            cur = np.where(a, cur + 1, 0)
            best = np.maximum(best, cur)
            nrun += a & ~prev
            prev = a
        out["kd_run_max"][rows] = best
        out["kd_run_n"][rows] = nrun
        out["kd_run_frac"][rows] = above.sum(1) / (int((~anc).sum()) or 1)

        ok = valid[:, :-1] & valid[:, 1:] & (~anc)[None, :-1] & (~anc)[None, 1:]
        if ok.any():
            code2 = safe[:, :-1] * 20 + safe[:, 1:]
            pair_code.append(code2[ok])
            pair_row.append(np.repeat(rows, L - 1)[ok.ravel()])

    c_all = c_anc + c_tcr
    out["pc1"] = c_all @ BASIS["pc1"]
    out["pc2"] = c_all @ BASIS["pc2"]
    for role, C in (("anchor", c_anc), ("tcr", c_tcr)):
        out[f"pc1_{role}"] = C @ BASIS["pc1"]
        out[f"pc2_{role}"] = C @ BASIS["pc2"]
        out[f"kf4_{role}"] = C @ BASIS["kf4"]
        out[f"mj_{role}"] = C @ BASIS["mj"]
    ntcr = np.maximum(c_tcr.sum(1), 1.0)
    out["para_tcr"] = (c_tcr @ BASIS["para"]) / ntcr
    out["para_sd_tcr"] = (c_tcr @ BASIS["para_sd"]) / ntcr
    cat = np.concatenate
    empty = np.empty(0, dtype=np.int64)
    return out, {"anchor": c_anc, "tcr": c_tcr, "n": n,
                 "pair_code": cat(pair_code) if pair_code else empty,
                 "pair_row": cat(pair_row) if pair_row else empty}


def apply_log_odds(counts: dict, source: str, weights) -> np.ndarray:
    """A fitted log-odds table applied to one count structure -- one value per peptide.

    ``anchor``/``tcr`` are dense (n, 20) and this is a matrix product. ``pair`` is the sparse pair
    list, and summing ``weights[code]`` per row with :func:`numpy.bincount` costs O(pairs) rather
    than O(n x 400) -- the difference between seconds and a gigabyte of temporaries on a corpus."""
    w = np.asarray(weights, dtype=np.float64)
    if source == "pair":
        return np.bincount(counts["pair_row"], weights=w[counts["pair_code"]],
                           minlength=counts["n"])[:counts["n"]]
    return counts[source] @ w


def design(peptides) -> np.ndarray:
    """The (n, k) design matrix in :func:`feature_names` order, before standardization."""
    feats, counts = encode(peptides)
    cols = [apply_log_odds(counts, FITTED[c], PARAMS["log_odds"][c]) if c in FITTED else feats[c]
            for c in feature_names()]
    return np.column_stack(cols)


def features(peptide: str) -> dict[str, float]:
    """The feature vector for one peptide, as a name -> value dict. For inspection; use
    :func:`score` (which takes a list) for anything at scale."""
    return dict(zip(feature_names(), design([peptide])[0].tolist()))


def score(peptides) -> np.ndarray:
    """Log-odds of immunogenic vs not, one per peptide. **Carries no prior.**

    Larger is more immunogenic. The training corpus's own base rate is divided out, so this is
    directly comparable across settings and composes with any prevalence via :func:`posterior`.

    Accepts a single string as well, and still returns an array of length 1 -- so a caller never has
    to branch on whether they passed one peptide or a million."""
    if isinstance(peptides, str):
        peptides = [peptides]
    peptides = list(peptides)
    if not peptides:
        return np.empty(0)
    X = design(peptides)
    st = PARAMS["standardizer"]
    Z = (X - np.asarray(st["mean"])) / np.asarray(st["std"])
    lo = PARAMS["logistic"]
    prev = PARAMS["prevalence"]
    return Z @ np.asarray(lo["coef"]) + lo["intercept"] - math.log(prev / (1.0 - prev))


def posterior(peptides, prior: float) -> np.ndarray:
    """``P(immunogenic | peptide)`` at an explicit ``prior``. Exact, because :func:`score` has none.

    ``prior`` has no default on purpose. The training corpus runs at ~3.2% positives, a viral
    proteome scan nearer 3.0e-3 and the NCI screen 4.2e-4 -- a default would silently pick one
    setting's base rate for every caller and overstate the rest by up to 75x."""
    if not 0.0 < prior < 1.0:
        raise ValueError(f"prior must be in (0, 1), got {prior!r}")
    z = score(peptides) + math.log(prior / (1.0 - prior))
    return 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))


def demo() -> None:
    """Self-check: run with ``python -m mhcmatch.complement``."""
    from . import posbayes

    # The role split is posbayes', residue for residue -- so the `aa` block is that model, and the
    # shipped table applied to these counts reproduces its score exactly.
    assert posbayes.AA == AA and tuple(posbayes.ANCHORS) == ANCHORS
    peps = ["GILGFVFTL", "SIINFEKL", "NLVPMVATV", "KRWIILGLNK", "RAKFKQLLA"]
    f, ct = encode(peps)
    for p, ca, cc in zip(peps, ct["anchor"], ct["tcr"]):
        r = posbayes.roles(len(p))
        wa, wt = np.zeros(20), np.zeros(20)
        for i, ch in enumerate(p):
            (wa if r[i] else wt)[AA.index(ch)] += 1
        assert (ca == wa).all() and (cc == wt).all(), p
    t = posbayes.table("human")
    mine = ct["anchor"] @ np.array(t["anchor"]) + ct["tcr"] @ np.array(t["tcrface"])
    assert max(abs(m - posbayes.llr(p)) for m, p in zip(mine, peps)) < 1e-9

    # The two role matrices partition the peptide, and whole-peptide sums are their total.
    assert (ct["anchor"].sum(1) + ct["tcr"].sum(1) == f["length"]).all()
    assert abs(f["pc1"][0] - (f["pc1_anchor"][0] + f["pc1_tcr"][0])) < 1e-9

    # The sparse pair list agrees with the dense matrix it replaces, row for row. The five anchors
    # are contiguous at each end, so an L-mer has L-5 TCR-facing positions in one block and exactly
    # L-6 adjacent pairs -- 3 for a 9-mer, 4 for a 10-mer.
    dense = np.zeros((len(peps), 400))
    np.add.at(dense, (ct["pair_row"], ct["pair_code"]), 1)
    assert dense.sum(1).tolist() == [len(p) - 6 for p in peps], dense.sum(1)
    wv = np.arange(400, dtype=float) / 400.0
    assert np.allclose(apply_log_odds(ct, "pair", wv), dense @ wv)

    # Composition features are permutation-invariant; runs are not. TCR-facing positions of a
    # 9-mer are 3..6, so IIDD is one run of 2 and IDID is two runs of 1 -- same composition.
    g, _ = encode(["AAAIIDDAA", "AAAIDIDAA"])
    assert (g["kd_run_max"][0], g["kd_run_n"][0]) == (2.0, 1.0)
    assert (g["kd_run_max"][1], g["kd_run_n"][1]) == (1.0, 2.0)
    assert g["kd_run_frac"][0] == g["kd_run_frac"][1]

    # Case, whitespace and non-standard residues.
    assert abs(score(["GILGFVFTL"])[0] - score([" gilgfvftl "])[0]) < 1e-12
    assert np.isfinite(score(["GILGFVFTX"])[0])
    assert len(score([])) == 0 and len(score("GILGFVFTL")) == 1

    # Batching must not change a score, and must be the reason to batch.
    a = score(peps)
    b = np.concatenate([score(peps[:2]), score(peps[2:])])
    assert np.allclose(a, b)

    # A lower prior lowers every posterior without touching the ranking.
    hi, lo = posterior(peps, PARAMS["prevalence"]), posterior(peps, 3.0e-3)
    assert (hi > lo).all() and (lo > 0).all()
    assert list(np.argsort(hi)) == list(np.argsort(lo)) == list(np.argsort(a))

    # The physics claim as a check: at fixed length, more hydrophobic TCR-facing residues score
    # higher -- the direction Chowell 2015 reports.
    assert score(["AAAIIIIAA"])[0] > score(["AAADDDDAA"])[0]

    print(f"ok - {len(feature_names())} features over {len(BLOCKS)} blocks, "
          f"fitted on {PARAMS['n']:,} rows ({PARAMS['arm']}); "
          f"score(GILGFVFTL) = {score(['GILGFVFTL'])[0]:+.4f}, "
          f"P@corpus = {posterior(['GILGFVFTL'], PARAMS['prevalence'])[0]:.4f}")


if __name__ == "__main__":
    demo()
