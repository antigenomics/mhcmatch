"""Recognition score for MHC epitopes: how immunogenic a peptide looks, given how it is presented.

Three heads, each fitted alone so that their fit criteria are comparable and each score is readable
on its own terms. The default is whichever wins BIC on the training arm, which is currently
``posbayes`` for both species:

``posbayes``
    Naive Bayes over amino-acid identity conditioned on **face** -- MHC-facing or TCR-facing --
    scored as a summed log-likelihood ratio. Two 20-cell tables, three parameters. Conditioning on
    the face rather than on absolute position is what keeps it well defined when peptide length
    varies, and it is why the two tables can disagree in sign without anything being told to flip.
    Pure numpy, no optional dependency, and the whole model prints in forty numbers.

``physchem_glm``
    Raw sums of the Kidera factors over each face. Length is not a separate feature: summing a
    constant-1 factor over a face gives that face's size, so it enters as :math:`KF_0` and the two
    face sizes add to the peptide length. 22 features, all interpretable.

``esm64_glm``
    64 principal components of a whole-peptide ESM2 pool. The most accurate head on mouse and the
    least explainable; needs ``pip install 'mhcmatch[esm]'``.

All three take what the caller knows, at whatever resolution they know it: peptide with both masks
given, peptide with an MHC (masks from the allele's layout), or peptide alone (the class-I default).

Coefficients come from ``chowell_iedb_full_matched`` -- the rebuilt Chowell corpus with negatives
resampled so the allele group carries no signal about the label, so no coefficient can be paid for
recognising which allele happened to be typed. Selection, transfer and the full metric matrix are in
``bench/results/shipped_models.md``.
"""

from __future__ import annotations

import functools
import json
import os

import numpy as np

__all__ = ["HEADS", "FITTED_HEADS", "SPECIES", "table", "default_head", "lowest_bic_head", "feature_names", "roles_for", "design",
           "score", "posterior", "embed", "mhc2_core", "score_mhc2", "MHC2_ANCHORS",
           "log_odds_table"]

SPECIES = ("human", "mouse")
#: ``complement`` is the six-block, 30-feature model of :mod:`mhcmatch.complement` and is the
#: **default**. The other three are fitted alone on one arm so their BIC is comparable to each
#: other; ``posbayes`` wins that comparison at three parameters. That is a statement about
#: parsimony on one corpus, not about which term to score with: in the integrated neoantigen fit
#: it is the six-block form that carries the recognition signal, and a 3-parameter head is a
#: different claim that must not be substituted for it silently. BIC across the three is still
#: reported by :func:`table`.
HEADS = ("complement", "posbayes", "physchem_glm", "esm64_glm")

#: The three heads with their own vendored artifacts. ``complement`` is served by its own module.
FITTED_HEADS = ("posbayes", "physchem_glm", "esm64_glm")
ESM_MODEL = "facebook/esm2_t33_650M_UR50D"
KIDERA = tuple(f"KF{i}" for i in range(1, 11))
#: 0-based P1/P4/P6/P9 within the register-anchored 9-mer core.
MHC2_ANCHORS = (0, 3, 5, 8)


def _data(name: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name)


@functools.lru_cache(maxsize=None)
def table(head: str = None, species: str = "human") -> dict:
    """Fitted parameters for one head. ``head=None`` resolves to the BIC-chosen default."""
    if species not in SPECIES:
        raise KeyError(f"no recognition table for {species!r}; have {', '.join(SPECIES)}")
    head = head or default_head(species)
    if head == "complement":
        raise KeyError("the 'complement' head has no vendored recognition table of its own; "
                       "it is mhcmatch.complement -- use complement.parameters(species)")
    if head not in HEADS:
        raise KeyError(f"unknown head {head!r}; have {', '.join(HEADS)}")
    with open(_data(f"recognition_{head}_mhc1_{species}.json")) as fh:
        p = json.load(fh)
    k = len(p["features"])
    if len(p["coef"]) != k or len(p["standardizer"]["mean"]) != k:
        raise ValueError(f"recognition/{head}/{species}: coefficients do not match features")
    return p


@functools.lru_cache(maxsize=1)
def _defaults() -> dict:
    with open(_data("recognition_default.json")) as fh:
        return json.load(fh)


def default_head(species: str = "human") -> str:
    """What :func:`score` uses when no head is named: the six-block ``complement`` model.

    :func:`lowest_bic_head` still reports the parsimony winner among the three separately fitted
    heads. They answer different questions -- BIC asks which head buys its parameters on one
    training arm, and the default asks which recognition term to score with.
    """
    return "complement"


def lowest_bic_head(species: str = "human") -> str:
    """The lowest-BIC head among the three with their own fitted tables."""
    return _defaults()["default"][species]


def feature_names(head: str = None, species: str = "human") -> list[str]:
    return list(table(head, species)["features"])


def log_odds_table(species: str = "human") -> dict:
    """``posbayes``'s two 20-cell tables as ``{face: {residue: log-odds}}``. The whole model."""
    t = table("posbayes", species)
    aa = t["alphabet"]
    return {face: dict(zip(aa, vals)) for face, vals in t["log_odds"].items()}


# ------------------------------------------------------------------ faces

def roles_for(peptides, mhc=None, anchors=None, store=None, cls="MHCI"):
    """Per-peptide boolean masks, ``True`` where the residue faces the MHC.

    Resolution order is explicit ``anchors``, then the allele's layout via ``store``, then the
    class-I default. Everything downstream reads this, so it is the one thing worth getting right.
    """
    out = []
    default = tuple(table("posbayes")["anchors"])
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


# ------------------------------------------------------------------ the three designs

@functools.lru_cache(maxsize=1)
def _kidera_basis():
    from .data import aa_tables
    from . import complement as CM
    tab = aa_tables.DESCRIPTORS["KIDERA"]
    # KF0 is the constant 1: summed over a face it is that face's size, so length is not separate
    return np.vstack([np.ones(20), [[tab[k][a] for a in CM.AA] for k in KIDERA]])


def _codes(peptides):
    from . import complement as CM
    ix = {a: i for i, a in enumerate(CM.AA)}
    return [np.array([ix[c] for c in p]) for p in peptides]


def _posbayes(peptides, roles, t):
    X = np.zeros((len(peptides), 2))
    tab = {0: np.asarray(t["log_odds"]["tcr"]), 1: np.asarray(t["log_odds"]["anchor"])}
    for n, (code, m) in enumerate(zip(_codes(peptides), roles)):
        for j, c in enumerate(code):
            f = 1 if m[j] else 0
            X[n, f] += tab[f][c]
    return X


def _physchem(peptides, roles):
    vec = _kidera_basis()
    X = np.zeros((len(peptides), 22))
    for n, (code, m) in enumerate(zip(_codes(peptides), roles)):
        m = np.asarray(m)
        a_i, t_i = code[m], code[~m]
        for j in range(11):
            X[n, 2 * j] = vec[j][a_i].sum() if a_i.size else 0.0
            X[n, 2 * j + 1] = vec[j][t_i].sum() if t_i.size else 0.0
    return X


@functools.lru_cache(maxsize=2)
def _esm():
    try:
        import torch
        from transformers import AutoTokenizer, EsmModel
    except ImportError as e:                              # pragma: no cover - env dependent
        raise ImportError("the esm64_glm head needs ESM2 embeddings: "
                          "pip install 'mhcmatch[esm]'. The default head does not.") from e
    dev = "mps" if torch.backends.mps.is_available() else (
        "cuda" if torch.cuda.is_available() else "cpu")
    tok = AutoTokenizer.from_pretrained(ESM_MODEL)
    mdl = EsmModel.from_pretrained(ESM_MODEL).to(dev).eval()
    return torch, tok, mdl, dev


def embed(peptides, batch: int = 256) -> np.ndarray:
    """Whole-peptide mean-pooled ESM2 embeddings, ``(n, 1280)``.

    Batched **within** length groups: slicing batches out of a length-sorted list and keeping only
    what matches the batch's first length silently leaves the rest unembedded.
    """
    torch, tok, mdl, dev = _esm()
    out = np.zeros((len(peptides), mdl.config.hidden_size), dtype=np.float32)
    by_len: dict[int, list[int]] = {}
    for i, p in enumerate(peptides):
        by_len.setdefault(len(p), []).append(i)
    for L, group in sorted(by_len.items()):
        for s in range(0, len(group), batch):
            idx = group[s:s + batch]
            enc = tok([peptides[i] for i in idx], return_tensors="pt",
                      add_special_tokens=True, padding=True).to(dev)
            with torch.no_grad():
                h = mdl(**enc).last_hidden_state[:, 1:L + 1, :]
            out[idx] = h.mean(1).float().cpu().numpy()
    return out


def _esm64(peptides, t):
    pca = np.load(_data("recognition_esm_pca.npz"))
    n = t["esm"]["n_components"]
    if pca["all_components"].shape[0] < n:
        raise ValueError(f"recognition_esm_pca.npz carries "
                         f"{pca['all_components'].shape[0]} components, the head needs {n}")
    return (embed(peptides) - pca["all_mean"]) @ pca["all_components"][:n].T


def design(peptides, species: str = "human", head: str = None, *, mhc=None, anchors=None,
           store=None, cls: str = "MHCI", roles=None) -> np.ndarray:
    """The design matrix for one head, in :func:`feature_names` order, before standardization."""
    head = head or default_head(species)
    t = table(head, species)
    peptides = list(peptides)
    if roles is None:
        roles = roles_for(peptides, mhc=mhc, anchors=anchors, store=store, cls=cls)
    if head == "posbayes":
        return _posbayes(peptides, roles, t)
    if head == "physchem_glm":
        return _physchem(peptides, roles)
    return _esm64(peptides, t)


def score(peptides, species: str = "human", head: str = None, **kw) -> np.ndarray:
    """Recognition log-odds. Higher is more likely to be immunogenic.

    Not a probability -- pass it through :func:`posterior` with a prior for the set being scored.
    """
    head = head or default_head(species)
    if head == "complement":
        from . import complement as _c
        return _c.score(peptides, species=species)
    t = table(head, species)
    X = design(peptides, species, head, **kw)
    m = np.asarray(t["standardizer"]["mean"])
    s = np.asarray(t["standardizer"]["std"])
    return float(t["intercept"]) + ((X - m) / s) @ np.asarray(t["coef"])


def posterior(peptides, prior: float, species: str = "human", head: str = None, **kw):
    """``sigmoid(score + logit(prior))``. ``prior`` is the immunogenic rate expected of the set."""
    if not 0.0 < prior < 1.0:
        raise ValueError("prior must be a probability")
    z = score(peptides, species, head, **kw) + float(np.log(prior / (1.0 - prior)))
    return 1.0 / (1.0 + np.exp(-z))


# ------------------------------------------------------------------ class II, and what it is not

_WARNED = False


def mhc2_core(peptides, register_start=None):
    """The register-anchored 9-mer core of each class-II peptide, and where it starts."""
    from .store import _mhc2_register
    cores, starts = [], []
    for i, p in enumerate(peptides):
        s = (register_start[i] if isinstance(register_start, (list, tuple)) else register_start) \
            if register_start is not None else _mhc2_register(p)
        if s is None or not 0 <= s <= len(p) - 9:
            cores.append(None)
            starts.append(None)
        else:
            cores.append(p[s:s + 9])
            starts.append(int(s))
    return cores, starts


def score_mhc2(peptides, species: str = "human", head: str = None, *, register_start=None,
               warn: bool = True):
    """Class-II score on an **MHC-I-fitted model**. Read this before using the number.

    There is no fitted class-II recognition model here, and this is not one: no class-II
    immunogenicity corpus with a usable negative set exists to fit against. This applies the class-I
    coefficients to the class-II binding core, with the groove-facing positions redefined as
    P1/P4/P6/P9 of the register-anchored 9-mer.

    Scoring the **core** rather than the whole peptide keeps every feature in the range the model was
    fitted on -- nine residues, face sizes of four and five. Scoring a 15-mer directly would put the
    face sizes far outside it. But the coefficients were fitted where the groove-facing residues are
    the termini and the TCR-facing ones a contiguous middle, and in class II the two faces
    interleave, so a coefficient learned on one geometry is being read on another. The register is
    also a heuristic unless supplied, and an error in the frame moves every residue between faces.

    Rank class-II peptides against each other with it. Do not compare the values with class-I
    scores, do not read them as calibrated, and say which model produced any number you report.
    ``nan`` where no 9-mer core can be assigned.
    """
    global _WARNED
    if warn and not _WARNED:
        import warnings
        warnings.warn(
            "recognition.score_mhc2 applies MHC-I-fitted coefficients to a class-II core: there is "
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
    out[ok] = score(sub, species, head, roles=roles)
    return out
