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
``motif``     Contiguity of the **hydropathy stretch**: longest run, number of runs and
              above-threshold fraction of hydrophobic TCR-facing residues. A run of 3-4 is
              a different object from the same residues scattered, and **no sum can express
              the difference**. :data:`KD_THRESHOLD` defines "hydrophobic" and
              :func:`encode` defines what breaks a run.
``aa``        Residue **identity** as a log-odds per amino acid per role, and the block that
              knows the peptide's **geometry**: a pooled anchor/TCR pair whose sum is exactly
              :func:`mhcmatch.posbayes.llr`, one pair per length bin, and a position key.
              Same shape in both classes; the length bins are 8/9/10/11+ at class I and
              quartiles of 11-25 at class II, and the position key is relative thirds of the
              TCR face at class I and **register zones** at class II.
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

Parameters are vendored per species in ``mhcmatch/data/complement_mhc1_{human,mouse}.json`` and
never refitted at import. The two hosts are **never pooled**: different MHC, different thymic
repertoires, so one fit across them is fitting a mixture. ``score(peps, species="mouse")`` selects.
Provenance, the block-by-block cross-validation, the corpus-transfer matrix and the size-matched
cross-species transfer are in the benchmark repo (``bench/results/complementarity.md``).

**Both classes, and they are two fits rather than one with a parameter.** ``score(peps)`` is class
I: the role split is P1-P3, PΩ-1, PΩ at fixed peptide positions. ``score(peps, cls="mhc2")`` takes
the anchors from the P1/P4/P6/P9 core of the floating 9-mer register
(:func:`mhcmatch.store.anchor_indices`) and reads its own vendored tables, fitted on 603,781 human
and 50,258 mouse class-II peptides from the IEDB export. Passing a class-II ligand to the class-I
path labels the wrong residues as anchors and returns a confident, wrong number, so the class is an
argument and never inferred from the length.

Parameters per class and species in ``complement_{mhc1,mhc2}_{human,mouse}.json``. Validation is in
``bench/results/complementarity.md`` and ``complementarity_mhc2.md``.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from importlib import resources
from statistics import median

import numpy as np

from . import immuno, ipred
from .data import aa_tables

__all__ = ["AA", "ANCHORS", "PARATOPE", "PARATOPE_CONTACT", "PARAMS", "TABLES", "SPECIES", "CLASSES", "BLOCKS", "burial",
           "BLOCKS_MHC2", "FITTED", "FITTED_MHC2", "ZONES", "MHC2_ZONES", "LENGTH_BINS",
           "MHC2_LEN_EDGES", "mhc2_length_bin",
           "blocks", "fitted", "mhc2_anchors", "table",
           "encode", "apply_log_odds", "design", "feature_names", "features", "score", "kidera_design", "kidera_names",
           "posterior", "parameters"]

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


#: The same potential marginalised over the residues that actually **contact** peptide, rather than
#: over the whole loop -- ``f_contact(a) = sum_{clonotypes, i} P(contact | locus, i, L) * 1[cdr3_i = a]``,
#: geometry from the 370 structures and identity from the same 28M clonotypes.
#:
#: TCRen is a *directed contact* potential, so the flat loop composition behind :data:`PARATOPE` is
#: the wrong measure for it, and wrong in a biasing direction: only **35.4%** of TRB loop residues
#: ever contact a peptide, and the germline V-encoded head and J-encoded tail carrying most of the
#: flat mass are exactly the positions the structures put at ``P(contact) = 0.00`` (TRB, length 12:
#: positions 1-2 and 9-12). Conditioning on contact **is** the germline-flank trim; it needs no
#: separate rule. Composition is never taken from the crystals -- 370 complexes are a poor estimate
#: of which residues sit mid-CDR3 -- only the geometry is.
#:
#: It is a materially different vector, not a rescaling: Spearman against :data:`PARATOPE` is
#: **+0.7549**, 19 of 20 residues change rank, and the spread widens from 0.2440 to 0.3255. A
#: structure-free second route restricted to the N-D-N insert agrees at Spearman +0.7173.
#:
#: Generated by ``bench/immuno/paratope_contact_basis.py``; ``bench/results/paratope_contact.md``
#: has the composition and ``paratope_contact_basis.md`` this vector. Opt in with
#: ``score(..., paratope="contact")``; the default stays ``"loop"`` so no recorded number moves.
PARATOPE_CONTACT = {
    "G": (-0.0177, 0.2136), "Q": (0.0050, 0.2787), "M": (0.0131, 0.3750), "I": (0.0233, 0.4399),
    "S": (0.0310, 0.2542), "P": (0.0365, 0.2677), "Y": (0.0370, 0.2607), "F": (0.0375, 0.3135),
    "L": (0.0396, 0.1878), "V": (0.0596, 0.2969), "E": (0.0623, 0.3612), "N": (0.0678, 0.3208),
    "R": (0.0805, 0.4688), "W": (0.0887, 0.4643), "T": (0.0965, 0.3768), "C": (0.1266, 0.5443),
    "A": (0.1288, 0.4503), "K": (0.1341, 0.4223), "D": (0.1674, 0.6565), "H": (0.3078, 0.5235),
}

#: The two MHC classes, which are two different constructions and not one with a parameter.
CLASSES = ("mhc1", "mhc2")

#: Length bins for the class-I ``aa`` block: 8, 9, 10 and **11+**. Closed at both ends, which is a
#: shippability property and not only a variance one -- a model with one table per *observed* length
#: cannot score a 12-mer at all, and real inputs contain those.
LENGTH_BINS = (8, 9, 10, 11)
#: Relative thirds of the class-I TCR-facing face. The same cell means the same fraction along the
#: peptide at every length, which is what the contact profile already does for its per-position
#: weights.
TCR_THIRDS = ("tcr_n", "tcr_m", "tcr_c")
#: Register-relative zones for class II: the residues **before** the 9-mer core, the non-anchor
#: residues **inside** it, and the residues **after** it.
#:
#: This was built as the class-II *replacement* for :data:`LENGTH_BINS`, on the reasoning that a
#: floating core makes total length uninformative. **The measurement says join, not replace**
#: (``bench/results/complementarity_mhc2.md``): the zones earn +0.0029 human / +0.0034 mouse AUROC
#: over the pooled pair, total length earns +0.0070 / +0.0159, and carrying both beats either on
#: both hosts and both metrics. The reasoning was right about the core and wrong about the ligand --
#: an 18-mer and a 13-mer with the same core do present the same residues, but a class-II ligand's
#: length is the length of its **flanks**, which is its own covariate and not a register question.
MHC2_ZONES = ("nflank", "core", "cflank")
#: Class-II total-length quartile edges, from the corpus's own distribution (11-25, median 15).
#: :data:`LENGTH_BINS` cannot be reused: it clamps to 11 and would put every class-II ligand in one
#: bin.
MHC2_LEN_EDGES = (14, 16, 19)


def mhc2_length_bin(L) -> int:
    """Which class-II length quartile (0-3) a ligand of length ``L`` falls in. Closed both ends."""
    return int(sum(int(L) >= e for e in MHC2_LEN_EDGES))
#: Which zone matrices a class partitions its TCR-facing face into. The pooled ``tcr`` matrix is
#: their sum in both cases, so the split costs no extra pass and ``posbayes`` stays recoverable.
ZONES = {"mhc1": TCR_THIRDS, "mhc2": MHC2_ZONES}


def length_bin(L) -> int:
    """Which of :data:`LENGTH_BINS` a peptide of length ``L`` falls in."""
    return min(max(int(L), LENGTH_BINS[0]), LENGTH_BINS[-1])


#: Blocks shared by both classes, in the order the benchmark's cumulative ablation adds them. Only
#: ``aa`` differs between the classes, because only ``aa`` is keyed on peptide geometry.
_SHARED_BLOCKS = {
    "phys": ["pc1", "pc2", "length"],
    "role": ["pc1_anchor", "pc2_anchor", "pc1_tcr", "pc2_tcr", "kf4_anchor", "kf4_tcr"],
    "pot": ["mj_anchor", "mj_tcr", "para_tcr", "para_sd_tcr"],
    "motif": ["kd_run_max", "kd_run_n", "kd_run_frac"],
}
#: Class-I feature blocks. ``BLOCKS`` keeps its name and its contents: it is the class-I layout and
#: every vendored class-I artifact is written against it.
BLOCKS = {
    **_SHARED_BLOCKS,
    "aa": (["aa_anchor", "aa_tcr"]
           + [f"aa_{r}{b}" for b in LENGTH_BINS for r in ("anchor", "tcr")]
           + [f"aa_{t}" for t in TCR_THIRDS]),
    "kmer": ["kmer_llr"],
}
#: Class-II feature blocks. The ``aa`` block has the same **shape** as class I's -- the pooled role
#: pair, a position key, and a length key -- and differs only in what the position key is: relative
#: thirds of the TCR face at class I, register zones at class II. See :func:`encode`.
BLOCKS_MHC2 = {
    **_SHARED_BLOCKS,
    "aa": (["aa_anchor", "aa_tcr"]
           + [f"aa_{z}" for z in MHC2_ZONES]
           + [f"aa_{r}L{b}" for b in range(len(MHC2_LEN_EDGES) + 1)
              for r in ("anchor", "tcr")]),
    "kmer": ["kmer_llr"],
}
#: Columns computed from a fitted log-odds table rather than from the peptide alone:
#: ``feature -> which count matrix it weights``. ``<matrix>@<bin>`` is that matrix restricted to the
#: rows of one length bin, so a per-length table costs a mask rather than another count matrix.
FITTED = {"aa_anchor": "anchor", "aa_tcr": "tcr", "kmer_llr": "pair",
          **{f"aa_{r}{b}": f"{r}@{b}" for b in LENGTH_BINS for r in ("anchor", "tcr")},
          **{f"aa_{t}": t for t in TCR_THIRDS}}
#: The class-II equivalent. ``<matrix>@<bin>`` here indexes :func:`mhc2_length_bin`, not
#: :func:`length_bin` -- :func:`encode` writes the right one into ``counts["bin"]`` per class.
FITTED_MHC2 = {"aa_anchor": "anchor", "aa_tcr": "tcr", "kmer_llr": "pair",
               **{f"aa_{z}": z for z in MHC2_ZONES},
               **{f"aa_{r}L{b}": f"{r}@{b}" for b in range(len(MHC2_LEN_EDGES) + 1)
                  for r in ("anchor", "tcr")}}


def blocks(cls: str = "mhc1") -> dict:
    """The feature blocks of one class: :data:`BLOCKS` or :data:`BLOCKS_MHC2`."""
    _check_cls(cls)
    return BLOCKS if cls == "mhc1" else BLOCKS_MHC2


def fitted(cls: str = "mhc1") -> dict:
    """The fitted-column map of one class: :data:`FITTED` or :data:`FITTED_MHC2`."""
    _check_cls(cls)
    return FITTED if cls == "mhc1" else FITTED_MHC2


def _check_cls(cls: str) -> None:
    if cls not in CLASSES:
        raise ValueError(f"unknown cls {cls!r} (expected one of {CLASSES})")

_SRC = "complement_{cls}_{species}.json"
#: Species with a fitted table. Both are the ``chowell_rebuilt`` arm of their own host -- human
#: 464,161 rows / 14,712 immunogenic, mouse 47,140 / 5,154 -- never pooled, because the two hosts
#: have different MHC and different thymic repertoires and a fit across them is fitting a mixture.
SPECIES = ("human", "mouse")


def _load(species: str = "human", cls: str = "mhc1") -> dict:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r} (expected one of {SPECIES})")
    _check_cls(cls)
    src = _SRC.format(cls=cls, species=species)
    ref = resources.files("mhcmatch.data").joinpath(src)
    if not ref.is_file():
        raise FileNotFoundError(
            f"no fitted {cls} / {species} complementarity model: {src} is not vendored. "
            f"Rebuild it with bench/neoag/complement{'' if cls == 'mhc1' else '_mhc2'}.py in the "
            f"benchmark repo, or use cls='mhc1', which always ships.")
    with ref.open() as fh:
        p = json.load(fh)
    k = len(p["features"])
    if len(p["logistic"]["coef"]) != k:
        raise ValueError(f"{src}: {len(p['logistic']['coef'])} coefficients for {k} features")
    if len(p["standardizer"]["mean"]) != k or len(p["standardizer"]["std"]) != k:
        raise ValueError(f"{src}: standardizer does not cover {k} features")
    if not 0.0 < p["prevalence"] < 1.0:
        raise ValueError(f"{src}: prevalence {p['prevalence']!r} is not a base rate")
    for name, s in p["log_odds_source"].items():
        want = 400 if s == "pair" else 20
        if len(p["log_odds"][name]) != want:
            raise ValueError(f"{src}: {name} has {len(p['log_odds'][name])} cells, want {want}")
        base = s.split("@")[0]
        if base not in ("anchor", "tcr", "pair", *ZONES[cls]):
            raise ValueError(f"{src}: {name} weights unknown count matrix {s!r}")
    if list(p["features"]) != [c for b in blocks(cls).values() for c in b]:
        raise ValueError(f"{src}: feature list is not the {cls} block layout")
    return p


#: The frozen **class-I** models, one per species: standardizer, linear head, the fitted log-odds
#: tables, and the EM / supervised Gaussian parameters kept for comparison. Class-I artifacts are
#: loaded eagerly because every caller that predates the class split expects them present at import.
TABLES: dict = {s: _load(s, "mhc1") for s in SPECIES}
#: The human class-I table, for callers that predate the species split.
PARAMS: dict = TABLES["human"]
#: ``(cls, species) -> parameters``, filled on first use. Class II is loaded lazily so that a
#: missing class-II artifact is an error at the call that wanted it, naming the class, rather than
#: an ImportError on ``import mhcmatch``.
_TABLES2: dict = {}


def table(species: str = "human", cls: str = "mhc1") -> dict:
    """The fitted parameters for one ``(species, cls)``."""
    _check_cls(cls)
    if cls == "mhc1":
        t = TABLES.get(species)
        if t is None:
            raise ValueError(f"unknown species {species!r} (expected one of {SPECIES})")
        return t
    key = (cls, species)
    if key not in _TABLES2:
        _TABLES2[key] = _load(species, cls)
    return _TABLES2[key]


def _scale_vec(tab: dict) -> np.ndarray:
    return np.array([tab[a] for a in AA], dtype=float)


def _basis(paratope: str = "loop") -> dict[str, np.ndarray]:
    """The residue vectors the encoder projects onto. ``paratope`` selects which marginalisation of
    the TCRen potential the ``pot`` block reads -- :data:`PARATOPE` or :data:`PARATOPE_CONTACT`."""
    tab = {"loop": PARATOPE, "contact": PARATOPE_CONTACT}.get(paratope)
    if tab is None:
        raise ValueError(f"unknown paratope={paratope!r} (expected 'loop' or 'contact')")
    rs = ipred.residue_scores()
    return {
        "pc1": np.array([rs[a][0] for a in AA]),
        "pc2": np.array([rs[a][1] for a in AA]),
        "kf4": _scale_vec(aa_tables.DESCRIPTORS["KIDERA"]["KF4"]),
        "mj": _scale_vec(aa_tables.MJ_PARTITION),
        "para": np.array([tab[a][0] for a in AA]),
        "para_sd": np.array([tab[a][1] for a in AA]),
        "kd": _scale_vec(aa_tables.HYDROPHOBICITY["KyteDoolittle"]),
    }


BASIS = _basis()
#: ``paratope -> basis``. Only the two ``para*`` vectors differ between them, but they are built
#: whole so the encoder reads one dict and never has to know which knob produced it.
_BASES = {"loop": BASIS, "contact": _basis("contact")}
#: The "hydrophobic" cut for the ``motif`` run features: the median of the Kyte-Doolittle scale
#: itself over the 20 standard residues (**-0.85**), so it is a property of the scale rather than a
#: constant tuned on any corpus. It admits ``ACFGILMSTV`` and excludes the other ten. Same rule as
#: :func:`mhcmatch.immuno._aggregate`.
#:
#: Taken from the plain dict with the stdlib rather than from :data:`BASIS` with ``numpy``, because
#: sphinx mocks ``numpy`` at doc-build time and a module-level ``float(np.median(...))`` then raises
#: on a Mock -- the whole module fails to import and its page renders empty.
KD_THRESHOLD = float(median(aa_tables.HYDROPHOBICITY["KyteDoolittle"].values()))


def feature_names(species: str = "human", cls: str = "mhc1") -> list[str]:
    """Column order of the design matrix, matching the vendored coefficients."""
    return list(table(species, cls)["features"])


def parameters(species: str = "human", cls: str = "mhc1") -> dict:
    """The fitted model as a plain dict (a copy of the vendored file)."""
    return json.loads(json.dumps(table(species, cls)))


@lru_cache(maxsize=1 << 16)
def mhc2_anchors(peptide: str, register: int | None = None) -> tuple:
    """0-based P1/P4/P6/P9 of a class-II ligand's 9-mer core, memoised.

    Delegates to :func:`mhcmatch.store.anchor_indices`, so the register is the same one
    ``decompose``, the logos and the signatures use. ``register`` pins the frame -- pass
    :meth:`mhcmatch.diffusion.AnchorModel.best_register` to annotate with the frame the model
    actually scored with, rather than the allele-agnostic heuristic.

    The heuristic register is an argmax over up to seventeen offsets **per peptide** and is the only
    non-vectorised step in :func:`encode`, so it is cached: a corpus repeats its peptides across
    alleles and folds, and the register does not depend on either.
    """
    from .store import anchor_indices
    return anchor_indices(peptide, "mhc2", register)


def _layout(peps, cls: str, registers) -> dict:
    """``(length, anchor indices) -> row numbers``: rows that share a residue layout.

    For class I the anchors are a function of the length alone, so this is exactly a grouping by
    length. For class II the register moves *within* a length, so each ``(length, register)`` is its
    own group -- which is why the grouping key is the anchor tuple and not the length.
    """
    groups: dict = {}
    for i, p in enumerate(peps):
        L = len(p)
        if L < 3:
            continue
        if cls == "mhc1":
            key = (L, tuple(sorted(a % L for a in ANCHORS)))
        else:
            key = (L, mhc2_anchors(p, None if registers is None else registers[i]))
        groups.setdefault(key, []).append(i)
    return groups


def _zone_positions(L: int, anc: np.ndarray, anchor_idx: tuple, cls: str) -> list:
    """The TCR-facing positions of one layout, partitioned into that class's :data:`ZONES`."""
    tcr_pos = [j for j in range(L) if not anc[j]]
    if cls == "mhc1":
        # Relative thirds: the same cell means the same fraction along the peptide at every length.
        k = len(tcr_pos)
        return [tcr_pos[int(k * f):int(k * g)] for f, g in ((0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1))]
    # Register-relative: before the core, inside it, after it. A ligand too short to carry a core
    # has no anchors and no flanks, and every residue is filed under `core`.
    s = anchor_idx[0] if anchor_idx else 0
    return [[j for j in tcr_pos if j < s],
            [j for j in tcr_pos if s <= j < s + 9],
            [j for j in tcr_pos if j >= s + 9]]


#: The chemistry half of the recognition axis, as a scale name. ``"Rose"`` is the shipped choice --
#: the average fractional area a residue buries on folding -- selected out of 576 candidates by the
#: BIC change it produced *inside the general model*, not by standalone AUROC. Any other name here
#: is exploratory: it re-parameterises a published result and should be reported as a comparison,
#: never substituted silently.
PHYS_SCALE = "Rose"


def phys_scale(scale=PHYS_SCALE) -> dict:
    """Resolve a residue scale by name to ``{residue: value}``.

    Accepts a key of :data:`mhcmatch.data.aa_tables.HYDROPHOBICITY` (45 scales, including the
    shipped ``"Rose"``), a ``"FAMILY:COMPONENT"`` key into
    :data:`~mhcmatch.data.aa_tables.DESCRIPTORS` (e.g. ``"KIDERA:KF4"``), or a ready dict over
    :data:`AA`. Raises rather than guessing on an unknown name."""
    from .data import aa_tables as _T
    if isinstance(scale, dict):
        miss = set(AA) - set(scale)
        if miss:
            raise ValueError(f"scale is missing {sorted(miss)}")
        return scale
    if scale in _T.HYDROPHOBICITY:
        return _T.HYDROPHOBICITY[scale]
    if ":" in scale:
        fam, comp = scale.split(":", 1)
        if fam in _T.DESCRIPTORS and comp in _T.DESCRIPTORS[fam]:
            return _T.DESCRIPTORS[fam][comp]
    raise ValueError(
        f"unknown scale {scale!r}. Use a HYDROPHOBICITY key (e.g. 'Rose'), a "
        f"'FAMILY:COMPONENT' descriptor key (e.g. 'KIDERA:KF4'), or a dict over the 20 residues.")


def burial(peptides, cls: str = "mhc1", scale=PHYS_SCALE, registers=None) -> list[float]:
    """``C_phys``: a residue scale summed over the **TCR-facing** positions.

    With the default ``scale="Rose"`` this is the shipped chemistry term -- the Rose burial
    propensity, i.e. the average fraction of a residue's surface buried on folding, summed over the
    face a receptor reads. It has **no fitted residue parameters**: the basis is imported, which is
    why it cannot memorise the corpus's cysteine gradient (correlation with per-peptide cysteine
    count +0.108, against +0.688 for the full fitted :func:`score`).

    ``scale=`` is for **exploration, not for scoring**. Passing another basis re-parameterises a
    result that was selected by BIC inside the general model over 576 candidates, so a number
    produced with a different scale is a comparison and must be reported as one. The obvious
    comparisons -- Kidera factors, against which the Chowell-family literature is usually
    written -- are reachable as ``"KIDERA:KF4"``, ``"KIDERA:KF2"`` and so on; on the neoantigen
    corpus they lose to ``"Rose"``.

    >>> burial(["GILGFVFTL"])[0] > 0
    True
    >>> abs(burial(["GILGFVFTL"], scale="KIDERA:KF4")[0]) >= 0
    True
    """
    import numpy as np
    tab = phys_scale(scale)
    vec = np.array([tab[a] for a in AA], dtype=float)
    _, counts = encode(list(peptides), cls, registers=registers)
    return (counts["tcr"].astype(float) @ vec).tolist()


def encode(peptides, cls: str = "mhc1", registers=None, positions: str = "mask",
           paratope: str = "loop") -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """``(features, counts)`` for an iterable of peptides -- the vectorised feature builder.

    Everything additive is a matrix product against a residue vector, so the peptide-side feature
    set is two (n, 20) count matrices (one per role) times seven property vectors. Run statistics
    need position order and use a loop over **positions** (at most 25 iterations), never over
    peptides. Peptides are grouped by residue layout so each group is one batch of array operations.

    ``counts`` carries ``anchor`` and ``tcr`` (n, 20) -- the matrices the ``aa`` log-odds tables
    weight -- plus that class's :data:`ZONES` matrices and the adjacent TCR-facing residue pairs as
    a **sparse pair list**, ``pair_code`` (which of 400 pairs) and ``pair_row`` (which peptide). A
    9-mer has at most 3 such pairs, so a dense (n, 400) matrix would be 99% zeros and would cost
    1.5 GB of temporaries per pass on a 500k-peptide corpus; :func:`apply_log_odds` sums the sparse
    form instead.

    **The two classes differ in one place and it is deliberate.** Class I is anchored at fixed
    peptide positions (:data:`ANCHORS`), so its geometry is a function of the peptide's *length*:
    the ``aa`` block carries one table per :data:`LENGTH_BINS` bin, and the TCR face is split into
    relative thirds. A class-II ligand is anchored by a 9-mer core that **floats** inside an
    11--25-mer, so its total length is the length of its flanking regions and carries no register
    information -- an 18-mer and a 13-mer with the same core present the same residues to the
    receptor. Binning a class-II table on total length would therefore split it on a variable that
    says nothing about the object being modelled. The class-II construction bins on the register
    instead (:data:`MHC2_ZONES`): before the core, inside it, after it. That both classes reach the
    same ``aa_anchor`` / ``aa_tcr`` pooled pair on top of it is what keeps
    :func:`mhcmatch.posbayes.llr` a strict special case of either.

    ``registers`` (class II only) pins each peptide's core to an explicit frame, element for element
    with ``peptides``; ``None`` uses the allele-agnostic heuristic register.

    **The hydropathy stretch, exactly.** A residue enters the ``motif`` block when all three hold:
    it is **TCR-facing** (an anchor is buried in the groove and is not part of any stretch the
    receptor reads), it is one of the 20 standard residues, and its **Kyte-Doolittle** value exceeds
    :data:`KD_THRESHOLD`. ``kd_run_max`` is then the longest run of consecutive such positions,
    ``kd_run_n`` the number of runs (counted as rising edges, so ``IIDI`` is two runs and ``IIDD``
    is one), and ``kd_run_frac`` the above-threshold count divided by the number of TCR-facing
    positions. An anchor **breaks** a run rather than bridging it: two stretches either side of a
    buried residue are two stretches.

    Non-standard residues (``X`` masks, ``B``/``J``/``O``/``U``/``Z``) contribute to no count matrix
    and to no property sum, but they still count toward ``length`` -- and in the ``motif`` block they
    behave exactly like a below-threshold residue rather than like a gap. A mask **breaks a run**
    (``AAAIIXIAA`` gives ``kd_run_max = 2``, the same as ``AAAIIDIAA``, not the 3 it would give if
    the mask were transparent) and it sits in ``kd_run_frac``'s denominator while never entering its
    numerator. That is the conservative reading of an unknown residue -- it is not evidence of a
    continued hydrophobic stretch -- and it is stated here because the opposite was once claimed.
    """
    _check_cls(cls)
    peps = [str(p).strip().upper() for p in peptides]
    n = len(peps)
    zone_names = ZONES[cls]
    c_anc = np.zeros((n, 20))
    c_zone = {z: np.zeros((n, 20)) for z in zone_names}
    # The `<matrix>@<bin>` machinery is shared; what `bin` MEANS is per class. Class I bins on the
    # peptide length clamped to 8-11; class II on the quartile of the 11-25 ligand range, because
    # LENGTH_BINS would put every class-II ligand in one bin.
    _bin = length_bin if cls == "mhc1" else mhc2_length_bin
    bins = np.array([_bin(len(p)) for p in peps], dtype=np.int16)
    if positions not in ("mask", "profile"):
        raise ValueError(f"unknown positions={positions!r} (expected 'mask' or 'profile')")
    if paratope not in _BASES:
        raise ValueError(f"unknown paratope={paratope!r} (expected 'loop' or 'contact')")
    basis = _BASES[paratope]
    # The TCR face read over the *positional contact profile* rather than the binary anchor mask.
    # Only the TCR face: the profile is a TCR<->peptide contact frequency and says nothing about
    # burial in the groove, so re-weighting the anchor face with it would be meaningless. The
    # binary `c_anc` / `c_tcr` are built regardless and are what the `aa` and `kmer` blocks read,
    # so those blocks are bit-identical under either setting.
    c_tcr_w = np.zeros((n, 20)) if positions == "profile" else None
    w_of = immuno.contact_profile(cls) if positions == "profile" else None
    pair_code: list = []
    pair_row: list = []
    out = {k: np.zeros(n) for k in ("kd_run_max", "kd_run_n", "kd_run_frac")}
    out["length"] = np.array([len(p) for p in peps], dtype=float)

    for (L, anchor_idx), idx in _layout(peps, cls, registers).items():
        rows = np.array(idx)
        code = np.array([[_AAI.get(ch, -1) for ch in peps[i]] for i in rows], dtype=np.int16)
        anc = np.zeros(L, dtype=bool)
        anc[list(anchor_idx)] = True
        valid = code >= 0
        safe = np.where(valid, code, 0)

        # The TCR face is accumulated per zone; the pooled `tcr` matrix is their sum, so the split
        # costs no extra pass and the pooled model stays exactly recoverable from it.
        masks = [("anchor", anc, c_anc)]
        for name, cols in zip(zone_names, _zone_positions(L, anc, anchor_idx, cls)):
            m = np.zeros(L, dtype=bool)
            m[cols] = True
            masks.append((name, m, c_zone[name]))
        for _, mask, mat in masks:
            sel = valid & mask[None, :]
            if sel.any():
                np.add.at(mat, (np.repeat(rows, L)[sel.ravel()], safe[sel]), 1)

        if c_tcr_w is not None:
            # `contact_profile` already zeroes sub-threshold positions and rescales the survivors
            # to mean 1, so a kept position weighs 1 exactly as it does under the binary mask and
            # the two encodings are on one scale. Its zeros ARE its anchor call -- on a class-I
            # 9-mer they are P1/P2/P3/POmega, which is `ANCHORS` without POmega-1, and the profile
            # gets there from geometry without being told anchors exist.
            wv = np.asarray(w_of(L), dtype=float)
            selw = valid & (wv > 0)[None, :]
            if selw.any():
                np.add.at(c_tcr_w, (np.repeat(rows, L)[selw.ravel()], safe[selw]),
                          np.broadcast_to(wv, (len(rows), L))[selw])

        # Hydrophobic runs, TCR-facing only. A buried anchor BREAKS a run rather than bridging it:
        # from the receptor's point of view two stretches either side of a buried residue are two
        # stretches, not one.
        # Under `positions="profile"` the face is the profile's non-zero set, which on a class-I
        # 9-mer is `ANCHORS` without POmega-1 -- so POmega-1 stops breaking runs and can join one.
        # That is the substantive difference between the two encodings for this block.
        face_mask = (~anc) if c_tcr_w is None else (np.asarray(w_of(L), dtype=float) > 0)
        above = valid & face_mask[None, :] & (basis["kd"][safe] > KD_THRESHOLD)
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
        if c_tcr_w is None:
            out["kd_run_frac"][rows] = above.sum(1) / (int((~anc).sum()) or 1)
        else:
            # weighted fraction: the profile rescales kept positions to mean 1, so this reduces to
            # the unweighted fraction when the footprint is flat
            wv2 = np.asarray(w_of(L), dtype=float)
            out["kd_run_frac"][rows] = (above * wv2[None, :]).sum(1) / (wv2.sum() or 1.0)

        ok = valid[:, :-1] & valid[:, 1:] & (~anc)[None, :-1] & (~anc)[None, 1:]
        if ok.any():
            code2 = safe[:, :-1] * 20 + safe[:, 1:]
            pair_code.append(code2[ok])
            pair_row.append(np.repeat(rows, L - 1)[ok.ravel()])

    c_tcr = sum(c_zone[z] for z in zone_names)
    c_all = c_anc + c_tcr
    out["pc1"] = c_all @ basis["pc1"]
    out["pc2"] = c_all @ basis["pc2"]
    face = c_tcr if c_tcr_w is None else c_tcr_w
    for role, C in (("anchor", c_anc), ("tcr", face)):
        out[f"pc1_{role}"] = C @ basis["pc1"]
        out[f"pc2_{role}"] = C @ basis["pc2"]
        out[f"kf4_{role}"] = C @ basis["kf4"]
        out[f"mj_{role}"] = C @ basis["mj"]
    ntcr = np.maximum(face.sum(1), 1.0)
    out["para_tcr"] = (face @ basis["para"]) / ntcr
    out["para_sd_tcr"] = (face @ basis["para_sd"]) / ntcr
    cat = np.concatenate
    empty = np.empty(0, dtype=np.int64)
    return out, {"anchor": c_anc, "tcr": c_tcr, **c_zone, "bin": bins, "n": n,
                 "pair_code": cat(pair_code) if pair_code else empty,
                 "pair_row": cat(pair_row) if pair_row else empty}


def apply_log_odds(counts: dict, source: str, weights) -> np.ndarray:
    """A fitted log-odds table applied to one count structure -- one value per peptide.

    ``anchor`` / ``tcr`` / the three :data:`TCR_THIRDS` are dense (n, 20) and this is a matrix
    product. ``pair`` is the sparse pair list, and summing ``weights[code]`` per row with
    :func:`numpy.bincount` costs O(pairs) rather than O(n x 400) -- the difference between seconds
    and a gigabyte of temporaries on a corpus.

    ``"<matrix>@<bin>"`` is that matrix restricted to one :data:`LENGTH_BINS` bin. The restriction
    is applied to the **result**, not by materialising a masked copy: a per-length table is then a
    length-n mask rather than another (n, 20) matrix, which is what makes eight of them affordable.
    """
    w = np.asarray(weights, dtype=np.float64)
    if source == "pair":
        return np.bincount(counts["pair_row"], weights=w[counts["pair_code"]],
                           minlength=counts["n"])[:counts["n"]]
    if "@" in source:
        base, b = source.split("@")
        return np.where(counts["bin"] == int(b), counts[base] @ w, 0.0)
    return counts[source] @ w


#: Index of cysteine, the residue :func:`score` can be asked to go blind to -- see :func:`_mask_cys`.
_CYS = _AAI["C"]


def _mask_cys(log_odds: dict, source: dict) -> dict:
    """The fitted log-odds tables with every cysteine cell zeroed -- a copy, never in place.

    **Why this exists.** :mod:`mhcmatch.posbayes` zeroes cysteine in every shipped table *by
    construction* and says why: the corpus's positives are synthetic assayed peptides and its
    negatives are mass-spectrometry-eluted ligands, and free cysteine is systematically
    under-recovered by MS, so the residue marks the *platform* rather than the label. This module
    was fitted without that mask, and the artifact is large: the vendored human class-I tables give
    cysteine ``+1.62`` to ``+2.63`` nats, the largest cell in all thirteen, against a mean
    ``|cell|`` of 0.24--0.41. One V->C substitution in ``GILGFVFTL`` moves :func:`score` by
    ``+2.90`` -- 1.8x the entire score of the unsubstituted epitope -- where
    :func:`mhcmatch.posbayes.llr` moves by ``+0.06``.

    The mask is **off by default** because every recorded number for this model was produced
    without it; turning it on is a different model and is reported as one.

    ``pair`` is 20x20 flattened, so both the first and the second residue of a dipeptide have to be
    zeroed, not just one.
    """
    out = {}
    for name, w in log_odds.items():
        v = np.asarray(w, dtype=np.float64).copy()
        if source[name] == "pair":
            v = v.reshape(20, 20)
            v[_CYS, :] = 0.0
            v[:, _CYS] = 0.0
            v = v.reshape(-1)
        else:
            v[_CYS] = 0.0
        out[name] = v
    return out


def _block_mask(features, cls: str, names) -> np.ndarray:
    """Boolean column mask selecting the features of ``names`` out of the full design."""
    b = blocks(cls)
    if isinstance(names, str):
        names = (names,)
    names = tuple(names)
    if not names:
        raise ValueError(f"blocks=() selects nothing; pass a subset of {tuple(b)} or None")
    unknown = [n for n in names if n not in b]
    if unknown:
        raise ValueError(f"unknown feature block(s) {unknown} (expected {tuple(b)})")
    want = {c for n in names for c in b[n]}
    return np.array([f in want for f in features], dtype=bool)


def design(peptides, species: str = "human", cls: str = "mhc1", registers=None,
           mask_cys: bool = False, positions: str = "mask",
           paratope: str = "loop") -> np.ndarray:
    """The (n, k) design matrix in :func:`feature_names` order, before standardization."""
    t = table(species, cls)
    fit = fitted(cls)
    lo = _mask_cys(t["log_odds"], t["log_odds_source"]) if mask_cys else t["log_odds"]
    feats, counts = encode(peptides, cls, registers, positions, paratope)
    cols = [apply_log_odds(counts, fit[c], lo[c]) if c in fit else feats[c]
            for c in t["features"]]
    return np.column_stack(cols)


KIDERA = tuple(f"KF{i}" for i in range(1, 11))


def kidera_names() -> list[str]:
    """Column names of :func:`kidera_design`, in order."""
    return [f"kf{i}_{r}" for i in range(1, 11) for r in ("anchor", "tcr", "all")]


def kidera_design(peptides, anchors=None, roles=None) -> np.ndarray:
    """All ten Kidera factors, each summed over the anchors, the TCR face and the whole peptide.

    The fitted model uses one of the ten -- ``kf4_anchor`` and ``kf4_tcr``, the hydropathy axis --
    because that is the factor the role split was established on, not because the other nine were
    compared and rejected. This builds the full set so they can be. Label-free: it reads only the
    vendored Kidera table and :data:`ANCHORS`, so it needs no fitted artifact and no per-fold refit.

    >>> kidera_design(["GILGFVFTL"]).shape
    (1, 30)
    """
    from .data import aa_tables
    tab = aa_tables.DESCRIPTORS["KIDERA"]
    vec = np.array([[tab[k][a] for a in AA] for k in KIDERA])      # (10, 20)
    idx = {a: i for i, a in enumerate(AA)}
    out = np.zeros((len(peptides), 30))
    for n, p in enumerate(peptides):
        L = len(p)
        if roles is not None:
            anc = {j for j in range(L) if roles[n][j]}
        else:
            anc = {a % L for a in (ANCHORS if anchors is None else anchors)}
        code = np.array([idx[c] for c in p])
        a_i = np.array([code[j] for j in range(L) if j in anc], dtype=int)
        t_i = np.array([code[j] for j in range(L) if j not in anc], dtype=int)
        for j in range(10):
            v = vec[j]
            out[n, 3 * j] = v[a_i].sum() if a_i.size else 0.0
            out[n, 3 * j + 1] = v[t_i].sum() if t_i.size else 0.0
            out[n, 3 * j + 2] = v[code].sum()
    return out


def features(peptide: str, species: str = "human", cls: str = "mhc1") -> dict[str, float]:
    """The feature vector for one peptide, as a name -> value dict. For inspection; use
    :func:`score` (which takes a list) for anything at scale."""
    return dict(zip(feature_names(species, cls), design([peptide], species, cls)[0].tolist()))


def score(peptides, species: str = "human", cls: str = "mhc1", registers=None,
          blocks=None, mask_cys: bool = False, positions: str = "mask",
          paratope: str = "loop") -> np.ndarray:
    """Log-odds of immunogenic vs not, one per peptide. **Carries no prior.**

    Larger is more immunogenic. The training corpus's own base rate is divided out, so this is
    directly comparable across settings and composes with any prevalence via :func:`posterior`.

    Accepts a single string as well, and still returns an array of length 1 -- so a caller never has
    to branch on whether they passed one peptide or a million.

    ``blocks`` restricts the score to a subset of :data:`BLOCKS`. The head is linear over a fixed
    block partition, so this is an **exact partial sum of the shipped score**, not a refit: the
    parts add back to the whole, up to the constant offset that belongs to no block. With
    ``off = logistic.intercept - log(prev / (1 - prev))`` from :func:`table`,

        score(p, blocks=("phys", "role", "pot", "motif")) + score(p, blocks=("aa", "kmer")) + off
            == score(p)

    exactly. The offset is dropped from a block score deliberately -- an intercept is not
    attributable to any block, and a constant cannot change a ranking. **Whole blocks only**:
    ``aa_tcr`` is fitted at ``-0.8796`` while every per-length ``aa_tcr{8,9,10,11}`` and every
    ``aa_tcr_{n,m,c}`` is positive, so the pooled column is collinear with its own decompositions
    and the fit pushed it negative to compensate. Dropping columns *within* a block without
    refitting inverts the dominant term; dropping a whole block cannot.

    ``mask_cys`` scores with cysteine zeroed out of the fitted log-odds tables, as
    :mod:`mhcmatch.posbayes` does by construction. See :func:`_mask_cys` for the size of the
    artifact and why the default is off.

    ``positions="profile"`` reads the TCR-facing *chemistry* over the **positional contact
    profile** (:func:`mhcmatch.immuno.contact_profile`: per-position TCR-peptide contact frequency
    from 8,062 contacts over 370 crystals, keyed by class and length) instead of the binary anchor
    mask. It touches the ``role``, ``pot`` and ``motif`` blocks only -- ``aa`` and ``kmer`` are
    count matrices keyed on a role and a weighted count is a different object, so they come back
    bit-identical. The anchor face is untouched for the same reason the TCR face is re-weighted:
    the profile measures contact with a *receptor* and carries no information about burial in the
    groove.

    ``paratope="contact"`` reads the ``pot`` block's TCRen projection off
    :data:`PARATOPE_CONTACT` -- the potential marginalised over the receptor residues that actually
    contact peptide -- instead of :data:`PARATOPE`'s flat whole-loop composition. **Independent of
    ``positions``**: one changes which peptide positions are read, the other which receptor
    residues the potential was averaged over, and they compose.

    **The coefficients are unchanged under either, so both are the shipped head reading a different
    encoding, not a refit.** Whether either should be refit is the question a benchmark answers;
    the defaults stay ``"mask"`` and ``"loop"``.
    """
    if isinstance(peptides, str):
        peptides = [peptides]
    peptides = list(peptides)
    if not peptides:
        return np.empty(0)
    t = table(species, cls)
    X = design(peptides, species, cls, registers, mask_cys, positions, paratope)
    st = t["standardizer"]
    Z = (X - np.asarray(st["mean"])) / np.asarray(st["std"])
    lo = t["logistic"]
    coef = np.asarray(lo["coef"])
    if blocks is not None:
        keep = _block_mask(t["features"], cls, blocks)
        return Z[:, keep] @ coef[keep]
    prev = t["prevalence"]
    return Z @ coef + lo["intercept"] - math.log(prev / (1.0 - prev))


def posterior(peptides, prior: float, species: str = "human", cls: str = "mhc1") -> np.ndarray:
    """``P(immunogenic | peptide)`` at an explicit ``prior``. Exact, because :func:`score` has none.

    ``prior`` has no default on purpose. The training corpus runs at ~3.2% positives, a viral
    proteome scan nearer 3.0e-3 and the NCI screen 4.2e-4 -- a default would silently pick one
    setting's base rate for every caller and overstate the rest by up to 75x."""
    if not 0.0 < prior < 1.0:
        raise ValueError(f"prior must be in (0, 1), got {prior!r}")
    z = score(peptides, species, cls) + math.log(prior / (1.0 - prior))
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

    # Both species tables load, are the same shape, and are genuinely different fits.
    for sp in SPECIES:
        assert len(table(sp)["logistic"]["coef"]) == len(feature_names(sp))
    assert feature_names("human") == feature_names("mouse")
    assert not np.allclose(score(peps, "human"), score(peps, "mouse"))
    for bad in ("rat", "Human", ""):
        try:
            score(peps, bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"score accepted species {bad!r}")

    # A lower prior lowers every posterior without touching the ranking.
    hi, lo = posterior(peps, PARAMS["prevalence"]), posterior(peps, 3.0e-3)
    assert (hi > lo).all() and (lo > 0).all()
    assert list(np.argsort(hi)) == list(np.argsort(lo)) == list(np.argsort(a))

    # The physics claim as a check: at fixed length, more hydrophobic TCR-facing residues score
    # higher -- the direction Chowell 2015 reports.
    assert score(["AAAIIIIAA"])[0] > score(["AAADDDDAA"])[0]

    print(f"ok - {len(feature_names())} features over {len(BLOCKS)} blocks, "
          + " / ".join(f"{s}: {table(s)['n']:,} rows" for s in SPECIES) + "; "
          f"score(GILGFVFTL) = {score(['GILGFVFTL'])[0]:+.4f}, "
          f"P@corpus = {posterior(['GILGFVFTL'], PARAMS['prevalence'])[0]:.4f}")


if __name__ == "__main__":
    demo()
