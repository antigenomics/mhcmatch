"""Recognition score for MHC-I epitopes: what a T-cell repertoire is likely to see.

This is the definitive recognition head. It answers one question -- given a peptide, and as much as
is known about how it is presented, how immunogenic does it look -- and it takes that knowledge at
whatever resolution the caller has:

* **peptide, MHC and both masks.** Nothing is guessed. Use this when the interface is measured, or
  when a class-II register or a structure says which residues face the groove.
* **peptide and MHC.** The masks are derived from the allele, through the same anchor layout the
  presentation model uses.
* **peptide alone.** The class-I default split ``(0, 1, 2, -2, -1)`` is applied. Pass a
  :class:`~mhcmatch.store.Store` as ``store=`` to have the best-presenting allele chosen first and
  its own layout used instead.

The design is 105 features: 20 amino-acid counts, length, the ten Kidera factors summed separately
over the MHC-facing and TCR-facing residues, and 32 principal components each of ESM2 embeddings
mean-pooled over those same two faces. Fitted per species on the rebuilt Chowell corpus, never
pooled. Selection of both the feature set and the training arm is recorded in
``bench/results/recognition_model.md`` of the benchmark repository -- in particular, resampling the
negatives to match population HLA usage was measured and **loses** (-0.019 human, -0.058 mouse on
held-out published deposits), so the unmatched arm is what ships.

The ESM components need ``torch`` and ``transformers``; everything else is numpy. Install with
``pip install 'mhcmatch[esm]'``. :func:`score` raises a clear error rather than silently dropping
the block, because a model missing a third of its design is not the model that was validated.
"""

from __future__ import annotations

import functools
import json
import os

import numpy as np

__all__ = ["SPECIES", "PARAMS", "table", "feature_names", "roles_for", "design", "score",
           "posterior", "embed", "mhc2_core", "score_mhc2", "MHC2_ANCHORS"]

_SRC = "recognition_mhc1_{}.json"
_PCA = "recognition_esm_pca.npz"
SPECIES = ("human", "mouse")
ESM_MODEL = "facebook/esm2_t33_650M_UR50D"
#: 0-based P1/P4/P6/P9 within the register-anchored 9-mer core.
MHC2_ANCHORS = (0, 3, 5, 8)


def _data(name: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name)


def _load(species: str) -> dict:
    with open(_data(_SRC.format(species))) as fh:
        p = json.load(fh)
    k = len(p["features"])
    for key in ("coef", "blocks"):
        if len(p[key]) != k:
            raise ValueError(f"recognition/{species}: {key} has {len(p[key])} of {k} entries")
    for key in ("mean", "std"):
        if len(p["standardizer"][key]) != k:
            raise ValueError(f"recognition/{species}: standardizer {key} does not match features")
    return p


PARAMS = {s: _load(s) for s in SPECIES}


def table(species: str = "human") -> dict:
    if species not in PARAMS:
        raise KeyError(f"no recognition table for {species!r}; have {', '.join(SPECIES)}")
    return PARAMS[species]


def feature_names(species: str = "human") -> list[str]:
    return list(table(species)["features"])


def roles_for(peptides, mhc=None, anchors=None, store=None, cls="MHCI"):
    """Per-peptide boolean masks, ``True`` where the residue faces the MHC.

    Resolution order is explicit ``anchors``, then the allele's own layout via ``store``, then the
    class-I default. The returned masks are what :func:`design` splits on, so a caller that wants a
    measured interface simply passes it here.
    """
    out = []
    default = tuple(table("human")["anchors"])
    for i, p in enumerate(peptides):
        L = len(p)
        if anchors is not None:
            a = anchors[i] if isinstance(anchors[0], (list, tuple, set, np.ndarray)) else anchors
            idx = {int(x) % L for x in a}
        elif store is not None and mhc is not None:
            allele = mhc[i] if isinstance(mhc, (list, tuple)) else mhc
            idx = {int(x) % L for x in store.anchor_indices(p, cls, allele)}
        else:
            idx = {x % L for x in default}
        out.append([j in idx for j in range(L)])
    return out


@functools.lru_cache(maxsize=2)
def _esm():
    try:
        import torch
        from transformers import AutoTokenizer, EsmModel
    except ImportError as e:                                   # pragma: no cover - env dependent
        raise ImportError(
            "the recognition score needs ESM2 embeddings: pip install 'mhcmatch[esm]'") from e
    dev = "mps" if torch.backends.mps.is_available() else (
        "cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained(ESM_MODEL)
    mdl = EsmModel.from_pretrained(ESM_MODEL).to(dev).eval()
    return torch, tok, mdl, dev


def embed(peptides, roles, batch: int = 256) -> tuple[np.ndarray, np.ndarray]:
    """Mean-pooled ESM2 embeddings over the MHC-facing and the TCR-facing residues.

    Returns ``(anchor, tcr)``, each ``(n, 1280)``. Batched by length so no padding is needed and the
    mask is one slice per batch.
    """
    torch, tok, mdl, dev = _esm()
    n = len(peptides)
    dim = mdl.config.hidden_size
    A = np.zeros((n, dim), dtype=np.float32)
    T = np.zeros((n, dim), dtype=np.float32)
    # Group by length first, then batch inside each group. Batching a length-sorted list and
    # keeping only the entries matching the batch's first length leaves the rest unembedded.
    by_len: dict[int, list[int]] = {}
    for i, p in enumerate(peptides):
        by_len.setdefault(len(p), []).append(i)
    for L, group in sorted(by_len.items()):
        for start in range(0, len(group), batch):
            idx = group[start:start + batch]
            enc = tok([peptides[i] for i in idx], return_tensors="pt",
                      add_special_tokens=True, padding=True).to(dev)
            with torch.no_grad():
                h = mdl(**enc).last_hidden_state[:, 1:L + 1, :]     # strip BOS/EOS
            for row, i in enumerate(idx):
                anc = [j for j in range(L) if roles[i][j]]
                tcr = [j for j in range(L) if not roles[i][j]]
                if anc:
                    A[i] = h[row, anc, :].mean(0).float().cpu().numpy()
                if tcr:
                    T[i] = h[row, tcr, :].mean(0).float().cpu().numpy()
    return A, T


def design(peptides, species: str = "human", *, mhc=None, anchors=None, store=None,
           cls: str = "MHCI", roles=None) -> np.ndarray:
    """The ``(n, 105)`` design matrix in :func:`feature_names` order, before standardization."""
    from . import complement as CM
    t = table(species)
    peptides = list(peptides)
    if roles is None:
        roles = roles_for(peptides, mhc=mhc, anchors=anchors, store=store, cls=cls)

    aa = CM.AA
    ix = {a: i for i, a in enumerate(aa)}
    comp = np.zeros((len(peptides), len(aa)))
    for n, p in enumerate(peptides):
        for c in p:
            comp[n, ix[c]] += 1
    length = np.array([[len(p)] for p in peptides], float)

    kid = CM.kidera_design(peptides, roles=roles)
    keep = [i for i, c in enumerate(CM.kidera_names()) if not c.endswith("_all")]
    kid = kid[:, keep]

    A, T = embed(peptides, roles)
    pca = np.load(_data(_PCA))
    npc = t["esm"]["n_components"]
    ea = (A - pca["anchor_mean"]) @ pca["anchor_components"][:npc].T
    et = (T - pca["tcr_mean"]) @ pca["tcr_components"][:npc].T
    return np.column_stack([comp, length, kid, ea, et])


def score(peptides, species: str = "human", **kw) -> np.ndarray:
    """Recognition log-odds, on the scale the model was fitted on.

    Higher is more likely to be immunogenic. The value is **not** a probability -- pass it through
    :func:`posterior` with a prior appropriate to the set being scored.
    """
    t = table(species)
    X = design(peptides, species, **kw)
    m = np.asarray(t["standardizer"]["mean"])
    s = np.asarray(t["standardizer"]["std"])
    return ((X - m) / s) @ np.asarray(t["coef"])


def posterior(peptides, prior: float, species: str = "human", **kw) -> np.ndarray:
    """``sigmoid(score + logit(prior))``. ``prior`` is the immunogenic rate expected of the set."""
    if not 0.0 < prior < 1.0:
        raise ValueError("prior must be a probability")
    z = score(peptides, species, **kw) + float(np.log(prior / (1.0 - prior)))
    return 1.0 / (1.0 + np.exp(-z))


# ------------------------------------------------------------------ class II, and what it is not

_WARNED = False


def mhc2_core(peptides, register_start=None):
    """The register-anchored 9-mer core of each class-II peptide, and where it starts.

    Uses the one-pass heuristic register in :mod:`mhcmatch.store` unless ``register_start`` pins the
    frame -- pass :meth:`mhcmatch.diffusion.AnchorModel.best_register` to annotate with the same
    frame a per-allele model scored on. Peptides shorter than nine residues yield ``(None, None)``.
    """
    from .store import _mhc2_register
    cores, starts = [], []
    for i, p in enumerate(peptides):
        if register_start is not None:
            s = register_start[i] if isinstance(register_start, (list, tuple)) else register_start
        else:
            s = _mhc2_register(p)
        if s is None or not 0 <= s <= len(p) - 9:
            cores.append(None)
            starts.append(None)
        else:
            cores.append(p[s:s + 9])
            starts.append(int(s))
    return cores, starts


def score_mhc2(peptides, species: str = "human", *, register_start=None, warn: bool = True):
    """Recognition score for class-II peptides, on an **MHC-I-trained model**. Read this first.

    There is no fitted class-II recognition model here, and this is not one. There is no class-II
    immunogenicity corpus with a usable negative set to fit against: the negatives every arm in this
    work uses are eluted self ligands paired with T-cell-assay positives, and that construction has
    not been built for class II. What this function does is apply the class-I coefficients to the
    class-II binding core, with the groove-facing positions redefined as P1/P4/P6/P9 of the
    register-anchored 9-mer rather than the class-I P2/P|Omega| pattern.

    Two reasons that is worth something rather than nothing. The design is mostly interface geometry
    -- Kidera factors and ESM embeddings pooled over the groove-facing and TCR-facing residues -- and
    that split is defined for class II too. And scoring the **core** rather than the whole peptide
    keeps every feature in the range the model was fitted on: nine residues, composition summing to
    nine, length fixed. Scoring a 15-mer directly would put `length` five standard deviations outside
    the fitted range and scale every count with it.

    Two reasons to distrust the number anyway. The coefficients were fitted where the groove-facing
    residues are the termini and the TCR-facing ones are a contiguous middle; in class II the
    groove-facing positions are interior and the faces interleave, so a coefficient learned on one
    geometry is being read on another. And the class-II register is itself a heuristic here, so an
    error in the frame moves every residue between the two faces.

    Use it to rank class-II peptides against each other. Do not compare the values to class-I scores,
    do not read them as calibrated, and do not report them without saying which model produced them.

    Returns ``nan`` for any peptide with no assignable 9-mer core.
    """
    global _WARNED
    if warn and not _WARNED:
        import warnings
        warnings.warn(
            "recognition.score_mhc2 applies MHC-I-trained coefficients to a class-II core: there is "
            "no fitted class-II model and no class-II corpus to fit one on. Ranking only.",
            UserWarning, stacklevel=2)
        _WARNED = True
    peptides = list(peptides)
    cores, _ = mhc2_core(peptides, register_start)
    out = np.full(len(peptides), np.nan)
    ok = [i for i, c in enumerate(cores) if c is not None]
    if not ok:
        return out
    sub = [cores[i] for i in ok]
    roles = [[j in MHC2_ANCHORS for j in range(9)] for _ in sub]
    out[ok] = score(sub, species, roles=roles)
    return out

