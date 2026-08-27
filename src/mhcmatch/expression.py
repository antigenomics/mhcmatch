"""Reference expression, keyed by normal tissue (GTEx) and by tumour type (TCGA).

A neoantigen ranker needs two different expression questions answered, and they are not the same
number, so this module never merges them:

**Is the source gene transcribed in the tumour's lineage?** The ranking read. A candidate from a gene
that is silent in that tissue is not presented however well it binds. Melanoma is ``SKCM``.

**Is it also transcribed in normal tissue?** The safety read. A gene expressed everywhere is a
toxicity risk, not a target.

Two key types, joined differently, both carrying their own provenance:

* ``key_type="gene"`` x a **GTEx tissue** -- per-tissue median TPM over the v11 bulk release. This is
  what imputes a missing TPM and what answers the safety question.
* ``key_type="peptide"`` x a **TCGA ``cancer_type``** -- the median ``expressionEB++`` over every TCGA
  sample in which that exact peptide was called. Keyed on the peptide rather than the gene
  deliberately: the TCGA source table carries ``ensp`` and no ENSP->symbol map ships with it, so a
  gene-level join would be a guess. The peptide-level join is exact, and it answers a question the
  gene-level number cannot -- *has this exact neoantigen been seen expressed in this tumour type*.

The table is fetched from the public HF dataset :data:`mhcmatch.store.PMHC_REPO` and cached, like
every other reference here.

    >>> from mhcmatch import expression
    >>> expression.lookup("PMEL", tissue="Skin - Sun Exposed (Lower leg)")   # doctest: +SKIP
    {'median_tpm': 44.1, 'q25_tpm': 21.6, 'q75_tpm': 78.9, 'n': 605, 'source': 'gtex'}

**Missing is encoded, never dropped.** :func:`impute` returns the reference value *and* a flag saying
whether it was observed or imputed, so a caller can carry a missing-indicator column instead of
discarding the candidate -- which is the standing rule for every partially-covered covariate here.
"""
from __future__ import annotations

import functools
import logging
import gzip
import os

__all__ = ["REFERENCE_FILE", "REFERENCE_TOIL_FILE", "MATRIX_FILE", "FLOORS_FILE",
           "SYNONYMS_FILE", "COLUMNS", "TUMOR_TISSUE", "TUMOR_TISSUE_APPROXIMATE",
           "C_MIN", "C_MAX", "MIN_SHARED", "MIN_COVERAGE", "GAMMA_MIN", "GAMMA_MAX",
           "matched_tissues", "resolve_context",
           "fetch_reference", "fetch_matrix", "fetch_synonyms", "load", "lookup", "impute", "tissues",
           "tumor_types", "safety_profile",
           "context_floor", "gene_level", "batch_scale"]

#: Path inside the HF dataset repo.
REFERENCE_FILE = "expression/reference_expression.tsv.gz"

#: The companion table in which **GTEx and TCGA are the same unit**, ``TPM`` on every row, because
#: both cohorts went through one RSEM pipeline. :data:`REFERENCE_FILE` cannot answer a question that
#: spans the two -- its GTEx half is TPM, its TCGA half is RSEM normalised counts, and both are
#: written into ``median_tpm`` -- so anything that compares a tumour type with a tissue, or takes a
#: scale from one to divide into the other, reads this one instead. Gene-keyed throughout: 53 GTEx
#: tissues under ``source="toil_gtex"`` and 33 TCGA study codes under ``source="toil_tcga"``.
#:
#: **Nothing here parses it, and it is not in the bootstrap set.** It is 38.6 MB and the record; the
#: scoring path reads :data:`MATRIX_FILE`, :data:`FLOORS_FILE` and :data:`SYNONYMS_FILE`, which are
#: 6.6 MB between them and carry the same numbers. Fetch it with
#: ``fetch_reference(file=REFERENCE_TOIL_FILE)`` and read it with polars when an analysis wants the
#: rows themselves.
REFERENCE_TOIL_FILE = "expression/reference_expression_toil.parquet"

#: Columns of the reference table, in file order.
COLUMNS = ("key", "key_type", "source", "context", "median_tpm", "q25_tpm", "q75_tpm", "n")


def fetch_reference(path: str | None = None, file: str = REFERENCE_FILE) -> str:
    """Local path to a reference table, downloading it from HF on first use.

    ``file`` selects which one -- :data:`REFERENCE_FILE` by default, or
    :data:`REFERENCE_TOIL_FILE` for the single-pipeline GTEx/TCGA table.

    ``$MHCMATCH_EXPRESSION`` overrides, for offline and cluster runs, matching how
    ``$MHCMATCH_PMHC`` overrides the presentation panel. Pointed at a **directory** it resolves
    ``file``'s basename inside it, so one setting serves every deposit here.

    Pointed at a **file** it overrides :data:`REFERENCE_FILE` and that alone -- whatever the file is
    called, which is the long-standing contract. The other deposits resolve beside it and fall
    through to the download when they are not there. One path cannot stand in for four different
    files, and letting it try is how a caller ends up handing a gzipped TSV to ``np.load``."""
    if path:
        return path
    env = os.environ.get("MHCMATCH_EXPRESSION")
    if env:
        if not os.path.isfile(env):
            return os.path.join(env, os.path.basename(file))
        if file == REFERENCE_FILE:
            return env
        beside = os.path.join(os.path.dirname(env), os.path.basename(file))
        if os.path.isfile(beside):
            return beside
    from huggingface_hub import hf_hub_download

    from .store import PMHC_REPO
    return hf_hub_download(repo_id=PMHC_REPO, repo_type="dataset", filename=file)


@functools.lru_cache(maxsize=8)
def _resolve(path: str | None, _env: str | None) -> str:
    """:func:`fetch_reference`, memoized on ``(path, $MHCMATCH_EXPRESSION)``.

    With neither set, ``fetch_reference`` goes through ``hf_hub_download``, which is a cache lookup
    and costs ~0.5 s -- fine once, ruinous per call, and :func:`safety_profile` is called once per
    gene inside a loop. ``_env`` is in the key rather than read inside so that changing the
    environment re-resolves instead of returning the previous file."""
    return fetch_reference(path)


@functools.lru_cache(maxsize=4)
def load(path: str | None = None) -> dict:
    """``{(key_type, key, context): {median_tpm, q25_tpm, q75_tpm, n, source}}``, read once and cached.

    Read with the stdlib rather than a dataframe library: the package's runtime dependencies are
    numpy and huggingface_hub, and a dict keyed by the join tuple is what every caller here wants
    anyway."""
    out: dict = {}
    with gzip.open(fetch_reference(path), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {c: head.index(c) for c in COLUMNS}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(head):
                continue
            try:
                rec = {k: float(f[ix[k]]) for k in ("median_tpm", "q25_tpm", "q75_tpm")}
                rec["n"] = int(float(f[ix["n"]]))
            except ValueError:
                continue
            rec["source"] = f[ix["source"]]
            out[(f[ix["key_type"]], f[ix["key"]], f[ix["context"]])] = rec
    return out


def lookup(key: str, tissue: str | None = None, tumor: str | None = None,
           path: str | None = None) -> dict | None:
    """Reference expression for a gene in a normal ``tissue`` or a peptide in a ``tumor`` type.

    Exactly one of ``tissue`` / ``tumor`` is required -- they are different measurements in different
    units from different studies, and silently falling back from one to the other is how a tumour
    abundance ends up reported as a normal-tissue TPM."""
    if (tissue is None) == (tumor is None):
        raise ValueError("pass exactly one of tissue= (GTEx, gene-keyed) or tumor= (TCGA, "
                         "peptide-keyed)")
    tbl = load(path)
    if tissue is not None:
        return tbl.get(("gene", key, tissue))
    return tbl.get(("peptide", key, tumor))


def impute(key: str, observed: float | None = None, tissue: str | None = None,
           tumor: str | None = None, path: str | None = None) -> tuple[float | None, bool]:
    """``(value, was_imputed)`` -- the observed value if there is one, else the reference median.

    Never returns "drop this row". A candidate with no expression measurement still has a source
    gene whose typical expression in that tissue is known, and the flag lets a model carry a
    missing-indicator rather than losing the sample."""
    if observed is not None and observed == observed:
        return float(observed), False
    rec = lookup(key, tissue=tissue, tumor=tumor, path=path)
    return (rec["median_tpm"], True) if rec else (None, True)


def tissues(path: str | None = None) -> list[str]:
    """Every GTEx tissue in the reference table."""
    return sorted({c for (kt, _, c) in load(path) if kt == "gene"})


def tumor_types(path: str | None = None) -> list[str]:
    """Every TCGA ``cancer_type`` in the reference table (``SKCM`` is melanoma)."""
    return sorted({c for (kt, _, c) in load(path) if kt == "peptide"})


#: **Which normal tissue is a tumour type's matched normal**, so the safety read can be asked without
#: the caller having to know that melanoma pairs with skin.
#:
#: The two vocabularies are different and neither is clinical, which is worth being explicit about:
#:
#: * the keys are **TCGA study abbreviations** (NCI GDC), a research nomenclature. ``CRC`` is the one
#:   exception -- TCGA itself has ``COAD`` and ``READ`` separately, and the source table merged them.
#: * the values are **GTEx ``SMTSD``** tissue names, GTEx's own controlled vocabulary.
#: * neither is ICD-O-3, SNOMED CT or OncoTree. Nothing here maps to a clinical coding system, and a
#:   pipeline that needs one has to bring its own crosswalk.
#:
#: Curated by organ correspondence against the 53 tissues actually present in the reference table, so
#: every value resolves. Ordered best-match first. **``HNSC`` is the weak one and is marked**: GTEx
#: has no head-and-neck mucosa, so minor salivary gland and oesophageal mucosa are the nearest
#: epithelia rather than the matched normal. Where a tumour has more than one plausible normal, all
#: of them are listed rather than one being picked silently -- ``SKCM`` against sun-exposed and
#: sun-protected skin is a different safety question in each.
TUMOR_TISSUE: dict[str, tuple[str, ...]] = {
    "BLCA": ("Bladder",),
    "BRCA": ("Breast - Mammary Tissue",),
    "CESC": ("Cervix - Ectocervix", "Cervix - Endocervix"),
    "CRC":  ("Colon - Transverse", "Colon - Sigmoid"),
    "GBM":  ("Brain - Cortex", "Brain - Frontal Cortex (BA9)"),
    "HNSC": ("Minor Salivary Gland", "Esophagus - Mucosa"),      # approximate -- see above
    "KICH": ("Kidney - Cortex",),
    "KIRC": ("Kidney - Cortex",),
    "KIRP": ("Kidney - Cortex",),
    "LIHC": ("Liver",),
    "LUAD": ("Lung",),
    "LUSC": ("Lung",),
    "OV":   ("Ovary", "Fallopian Tube"),
    "PAAD": ("Pancreas",),
    "PRAD": ("Prostate",),
    "SKCM": ("Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"),
    "STAD": ("Stomach",),
    "THCA": ("Thyroid",),
    "UCEC": ("Uterus",),
}
#: Tumour types whose matched normal is an approximation rather than the same organ.
TUMOR_TISSUE_APPROXIMATE = ("HNSC",)


def matched_tissues(tumor: str) -> tuple[str, ...]:
    """The GTEx tissue(s) that are ``tumor``'s matched normal, best match first.

    ``()`` for a tumour type with no entry, which is the honest answer -- a wrong matched normal
    turns the safety read into a confident wrong one. See :data:`TUMOR_TISSUE` for what the two
    vocabularies are and what they are not."""
    return TUMOR_TISSUE.get(tumor.strip().upper(), ())


@functools.lru_cache(maxsize=4)
def _by_gene(resolved: str) -> dict:
    """``{gene: [(tissue, median_tpm)]}``, sorted descending -- one pass over the table.

    :func:`safety_profile` used to scan all 5,586,792 rows per call, **511 ms each**, and its callers
    ask per gene inside a loop: :func:`mhcmatch.mimicry.safety` once per mimic hit, and
    :func:`mhcmatch.vector.self_origin_risk` once per register of every candidate. A thousand-peptide
    screen therefore spent over half an hour in a linear scan it could pay once. Indexed: 0.1 us.

    **Keyed on the resolved file, not on the ``path`` argument.** Both are ``None`` when
    ``$MHCMATCH_EXPRESSION`` points somewhere else, so an argument-keyed cache would keep serving the
    previous table after the environment changed -- silently, and with plausible numbers."""
    out: dict = {}
    for (kt, key, ctx), row in load(resolved).items():
        if kt == "gene":
            out.setdefault(key, []).append((ctx, row["median_tpm"]))
    for v in out.values():
        v.sort(key=lambda x: -x[1])
    return out


#: A context with fewer gene rows than this does not estimate a percentile of its transcriptome.
_MIN_GENES = 30
_POOLED_SEEN: set = set()
_LOG = logging.getLogger(__name__)


@functools.lru_cache(maxsize=32)
def _tissue_quantile(path: str | None, tissues: tuple, q: float) -> tuple:
    """``(q-th percentile of non-zero gene medians in ``tissues``, n genes)``; ``nan`` if too few.

    Nearest-rank on the sorted values, so the result is exactly one observed abundance and two
    runs cannot disagree in the last bit. ``load`` is itself cached on ``path``, so this walks the
    table once per distinct ``(path, tissues, q)`` and is a dict lookup after that."""
    vals = []
    for (kt, _key, ctx), row in load(path).items():
        if kt == "gene" and (not tissues or ctx in tissues):
            v = row.get("median_tpm")
            if v is not None and v > 0:
                vals.append(v)
    if len(vals) < _MIN_GENES:
        return float("nan"), len(vals)
    vals.sort()
    i = min(len(vals) - 1, max(0, int(round(q * (len(vals) - 1)))))
    return float(vals[i]), len(vals)


def tissue_floor(tumor: str | None = None, tissue: str | None = None, q: float = 0.25,
                 path: str | None = None, detail: bool = False):
    """The abundance floor ``c`` for :func:`mhcmatch.rank.expr_level`, in TPM.

    The ``q``-th percentile of non-zero median abundance over **every gene** in a set of normal
    tissues -- the tumour's matched normal(s), or ``tissue`` named directly, or the whole reference
    pooled when neither is given.

    **This reads the gene half of the reference and never the peptide half.** The table holds two
    kinds of row and they answer different questions::

        gene rows      4,911,764   GTEx + HPA consensus, 104 normal tissues
                                   -> a whole transcriptome, so a percentile over it means something
        peptide rows   1,700,779   TCGA, 19 tumour types
                                   -> one peptide's abundance in one tumour type

    A floor is a statement about where a *transcriptome* stops resolving, so only the gene half can
    supply it. The peptide half supplies the numerator instead -- what a given candidate is
    transcribed at -- and the two must not be crossed.

    ``tumor`` is a TCGA code (``SKCM``, ``BLCA``, ``LUAD``); :func:`matched_tissues` maps it to its
    normal tissue(s). A code with no mapping falls back to the pooled reference, which is the honest
    answer when the origin is unknown. A ``tissue`` **named explicitly and not found is an error**,
    not a silent fallback: a typo would otherwise return a plausible number from the wrong
    distribution.

    With ``detail=True`` returns ``{"floor", "contexts", "n_genes", "pooled"}`` instead of the bare
    value, so a caller can see which tissues were used and whether it fell back.

    >>> round(tissue_floor(tumor="SKCM"), 4)                     # doctest: +SKIP
    0.1498
    >>> tissue_floor(tumor="SKCM", detail=True)["contexts"]      # doctest: +SKIP
    ('Skin - Sun Exposed (Lower leg)', 'Skin - Not Sun Exposed (Suprapubic)')
    """
    if not 0.0 < q < 1.0:
        raise ValueError(f"tissue_floor: q must be in (0, 1), got {q!r}")
    ts: tuple = ()
    if tissue:
        ts = (tissue,)
    elif tumor:
        ts = tuple(matched_tissues(tumor))
        if not ts:
            _log_pooled(tumor)
    v, n = _tissue_quantile(path, ts, q)
    if v != v and tissue:                                        # NaN-safe: named, and not found
        raise ValueError(
            f"tissue_floor: tissue {tissue!r} has fewer than {_MIN_GENES} gene rows in the "
            f"reference. Check the spelling against `mhcmatch expression --list-contexts`; "
            "falling back silently would return a number from the wrong distribution.")
    pooled = not ts or v != v
    if v != v:
        v, n = _tissue_quantile(path, (), q)
    if not detail:
        return v
    return {"floor": v, "contexts": ts, "n_genes": n, "pooled": pooled}


def _log_pooled(tumor: str) -> None:
    """A tumour code with no matched normal gets the pooled floor, and says so once."""
    if tumor not in _POOLED_SEEN:
        _POOLED_SEEN.add(tumor)
        _LOG.info("tissue_floor: %r has no matched normal tissue; using the pooled reference. "
                  "`mhcmatch expression --list-contexts` prints the tumour vocabulary.", tumor)


def safety_profile(gene: str, top: int = 10, path: str | None = None) -> list[tuple[str, float]]:
    """``[(tissue, median_tpm)]`` for a gene across normal tissues, highest first.

    The safety read: a target expressed only in the tumour's lineage is very different from one
    expressed in heart and lung too, and the ranking score alone does not show that."""
    return _by_gene(_resolve(path, os.environ.get("MHCMATCH_EXPRESSION"))).get(gene, [])[:top]


# ------------------------------------------------------------------ the single-pipeline reference

#: The whole single-pipeline table as a dense **gene x context float32 matrix**, which is what every
#: lookup on the scoring path reads. Measured against the same question asked of
#: :data:`REFERENCE_TOIL_FILE`, which builds a five-million-entry dictionary of dictionaries:
#: **0.05 s against 5.20 s to load, and 29 MB resident against 3,168 MB** -- 100x the speed on 1/109
#: of the memory. Both return the same floor to four decimals, which is how the matrix is checked.
MATRIX_FILE = "expression/toil_matrix.npz"

#: The same contexts' abundance floors at three quantiles, human-readable, 88 rows. Nothing on the
#: scoring path parses it -- :func:`context_floor` computes from the matrix, and this file is what a
#: table or a caption cites.
FLOORS_FILE = "expression/toil_floors.tsv"

#: **The floor is clamped to this range, in TPM.** Measured over all 86 contexts it runs 0.1000
#: (whole blood) to 0.4000 (testis), so the clamp is far wider than anything real and exists only so
#: that a degenerate input -- an empty context, a table read wrong, a caller passing a filter of zero
#: -- cannot produce a floor of 0 and a division that is not one.
C_MIN, C_MAX = 0.05, 2.0

#: A batch scale needs an *unconditioned* view of the transcriptome, and these two guards are what
#: enforce it. Absolute floor first: fewer shared genes than this is not an estimate.
MIN_SHARED = 1000
#: ...and then the one that actually matters -- the shared genes as a fraction of the genes the
#: reference context has switched on. A TCGA context expresses 23,508 to 30,561 genes (median
#: 26,518), a whole-transcriptome profile covers 0.6 to 0.9 of one, and **no candidate list can
#: reach 0.5**: the largest screen in the fitting corpus manages 4,772 shared genes, or 0.18.
#: Measured on the three screens carrying their own RNA-seq, a scale taken from their candidate
#: lists comes out at 1.78, 2.18 and 3.15 -- all above 1, none of them a unit difference. A
#: mutation reaches a candidate list only if it was seen in RNA, so a candidate's abundance is
#: conditioned on having been detected, and the ratio measures that conditioning instead of the
#: library. Counting more candidates cannot fix it, which is why the guard is coverage.
MIN_COVERAGE = 0.5
#: The magnitude a scale is allowed to reach. Wide on purpose: raw counts against TPM differ by
#: three orders of magnitude or more, and absorbing exactly that is what the estimate is for. These
#: bounds catch a non-finite or structurally broken result, **not** an unusual unit -- a clamp tight
#: enough to call 1000x implausible would reject the commonest case it exists to handle.
GAMMA_MIN, GAMMA_MAX = 1e-6, 1e6

#: **How a free-text origin becomes a context.** 232 rows, long format, one target per row, each
#: carrying where it came from: the TCGA study codes and GTEx tissues present in the matrix, Xena's
#: ``_primary_site`` and ``detailed_category`` joined to the study code on the sample barcode, and a
#: short list of curated spelling variants marked ``curated:spelling``. Derived rather than written
#: by hand because an organ maps to more than one study more often than not -- ``Lung`` is ``LUAD``
#: *and* ``LUSC``, ``Kidney`` is three -- and a crosswalk that picks one is silently wrong for the
#: rest. Kept as a deposit rather than a module constant so it can be corrected without a release.
SYNONYMS_FILE = "expression/context_synonyms.tsv"


def fetch_matrix(path: str | None = None) -> str:
    """Local path to :data:`MATRIX_FILE`, downloading it from HF on first use."""
    return fetch_reference(path, file=MATRIX_FILE)


def fetch_synonyms(path: str | None = None) -> str:
    """Local path to :data:`SYNONYMS_FILE`, downloading it from HF on first use."""
    return fetch_reference(path, file=SYNONYMS_FILE)


@functools.lru_cache(maxsize=2)
def _synonyms(path: str | None = None) -> dict:
    """``{lowercased alias: {"tcga_code": (...), "gtex_context": (...)}}``, read once and cached.

    An alias present under more than one ``alias_kind`` is merged: ``Lung`` is both a GTEx tissue
    and an organ, and an origin reading "Lung" means the lung, whichever table it was written
    against."""
    out: dict = {}
    with open(fetch_synonyms(path), encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ia, ik, it = head.index("alias"), head.index("target_kind"), head.index("target")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ia, ik, it):
                continue
            d = out.setdefault(f[ia].strip().lower(), {"tcga_code": [], "gtex_context": []})
            if f[ik] in d and f[it] not in d[f[ik]]:
                d[f[ik]].append(f[it])
    return {k: {kk: tuple(vv) for kk, vv in v.items()} for k, v in out.items()}


@functools.lru_cache(maxsize=2)
def _matrix(path: str | None = None) -> tuple:
    """``(gene index, context index, values, n_samples)``, read once and cached.

    Context keys are ``"<source>|<context>"``, so ``toil_tcga|SKCM`` and ``toil_gtex|Liver`` cannot
    be confused for one another by a caller holding a bare name."""
    import numpy as np

    # `allow_pickle` stays off. The deposit stores its label arrays as fixed-width unicode for
    # exactly this reason: a pickled array in a file fetched over the network executes on load.
    z = np.load(fetch_matrix(path))
    gi = {str(g): i for i, g in enumerate(z["genes"])}
    ci = {str(c): i for i, c in enumerate(z["contexts"])}
    return gi, ci, z["values"], z["n_samples"]


def _gtex_contexts(path: str | None = None) -> dict:
    """``{lowercased tissue name: the name as the matrix spells it}``.

    Spelling is the reason this exists rather than a set. Toil writes ``Skin - Sun Exposed (Lower
    Leg)`` and ``Brain - Frontal Cortex (Ba9)`` where GTEx v11 writes ``(Lower leg)`` and ``(BA9)``,
    so a tissue name carried from one table to the other matches on nothing at all."""
    _, ci, _, _ = _matrix(path)
    return {c.split("|", 1)[1].lower(): c.split("|", 1)[1]
            for c in ci if c.startswith("toil_gtex|")}


def resolve_context(text: str, path: str | None = None, approximate: bool = True,
                    detail: bool = False):
    """``(TCGA study codes, GTEx tissue names)`` for a free-text origin.

    Accepts what a submission actually carries: a study code (``SKCM``), an organ
    (``liver``, ``Lung``), or a GTEx tissue verbatim (``Skin - Sun Exposed (Lower Leg)``). Matching
    is case-insensitive throughout, and an organ that is more than one study returns all of them
    rather than one of them.

    **An unrecognised string raises.** A tumour type that silently became the pooled reference would
    return a plausible number computed from the wrong distribution, and nothing downstream could
    tell. An origin that is genuinely unknown is expressed by passing nothing, not by passing a word
    that does not resolve.

    >>> resolve_context("liver")                                  # doctest: +SKIP
    (('LIHC',), ('Liver',))
    >>> resolve_context("lung")                                   # doctest: +SKIP
    (('LUAD', 'LUSC'), ('Lung',))
    """
    s = str(text or "").strip()
    if not s:
        raise ValueError("resolve_context: empty origin. Pass None to mean 'unknown' -- an empty "
                         "string is a missing value that looks like a request.")
    syn = _synonyms(path)
    hit = syn.get(s.lower())
    exact = hit is not None
    if hit is None and approximate:
        hit, matched = _by_organ(s, syn)
        if hit is not None:
            _log_approximate(s, matched)
    if hit is None:
        raise ValueError(
            f"resolve_context: {s!r} is not a TCGA study code, an organ, a disease name or a GTEx "
            "tissue, and no organ name occurs inside it. Check it against "
            "`mhcmatch expression --list-contexts`; resolving it to the pooled reference would "
            "return a number from the wrong distribution with no way to tell it had happened.")
    if detail:
        return {"tcga_codes": hit["tcga_code"], "gtex_contexts": hit["gtex_context"],
                "exact": exact}
    return hit["tcga_code"], hit["gtex_context"]


def _by_organ(text: str, syn: dict) -> tuple:
    """The longest organ name occurring inside ``text``, and what it resolves to.

    A registry writes its own study names and they are not the ones any deposit here happens to
    carry: ``Kidney Renal Clear Cell Carcinoma`` is not Xena's ``Kidney Clear Cell Carcinoma`` and
    matches no alias exactly. The organ is still named in it, so it is read out rather than the row
    being lost -- the longest alias wins, so ``Colon`` does not pre-empt a longer phrase containing
    it, and the result is logged as approximate.

    This never invents a context. It only fires on an alias already in the table, and a string with
    no organ in it -- ``TIL``, ``PBMC``, ``healthy`` -- still raises, which is correct: those are not
    tumour types and resolving them to one would be worse than failing."""
    low = " " + " ".join(str(text).lower().replace("-", " ").split()) + " "
    best = None
    for alias in syn:
        if len(alias) < 4 or f" {alias} " not in low:
            continue
        if best is None or len(alias) > len(best):
            best = alias
    return (syn[best], best) if best else (None, None)


_APPROX_SEEN: set = set()


def _log_approximate(text: str, alias: str) -> None:
    """An origin resolved by organ rather than exactly says so once, not once per candidate."""
    if text not in _APPROX_SEEN:
        _APPROX_SEEN.add(text)
        _LOG.info("expression: %r matched no context exactly; resolved on the organ %r. "
                  "The floor and the matched normal are that organ's.", text, alias)


def _matched_toil(code: str, gt: dict) -> tuple:
    """:func:`matched_tissues` mapped onto the spellings this matrix actually uses."""
    return tuple(gt[t.lower()] for t in matched_tissues(code) if t.lower() in gt)


def _floor_from(keys: tuple, q: float, path: str | None = None) -> tuple:
    """``(q-th percentile of non-zero levels pooled over ``keys``, n genes)``; ``nan`` if too few."""
    import numpy as np

    _, ci, V, _ = _matrix(path)
    cols = [ci[k] for k in keys if k in ci]
    if not cols:
        return float("nan"), 0
    v = V[:, cols].reshape(-1)
    v = v[v > 0]
    if v.size < _MIN_GENES:
        return float("nan"), int(v.size)
    return float(np.quantile(v, q, method="nearest")), int(v.size)


def context_floor(tumor: str | None = None, tissue: str | None = None, q: float = 0.25,
                  path: str | None = None, clamp: bool = True, prefilter: float = 0.0,
                  detail: bool = False):
    """The abundance floor ``c`` for :func:`mhcmatch.rank.expr_level`, in TPM.

    The ``q``-th percentile of non-zero median abundance over **every gene** in a context, taken
    from the table where GTEx and TCGA share a pipeline, so a tumour type and a tissue are the same
    unit and the two can be mixed without a conversion.

    Resolution order, and it prefers the tumour deliberately::

        tumor  -> that TCGA study's own transcriptome        SKCM 0.1600, LUAD 0.2000 TPM
        tissue -> that GTEx tissue's transcriptome           Lung 0.3500, Liver 0.1800 TPM
        neither-> the pooled TCGA reference                  0.1800 TPM

    **A tumour's floor is not its matched normal's.** Measured across the pairs, a study sits at
    roughly half the tissue it arose in -- SKCM 0.1600 against skin 0.3050, BLCA 0.1700 against
    bladder 0.3600, LUAD 0.2000 against lung 0.3500. Scoring a tumour candidate against a normal
    floor puts the whole term about one unit low.

    ``tumor`` may be a study code or an organ; :func:`resolve_context` decides, and an organ that is
    several studies pools them. ``prefilter`` is an expression cut the candidates already passed, in
    TPM, and raises the floor to meet it -- a filter removes the range this term resolves. ``clamp``
    holds the result inside :data:`C_MIN`--:data:`C_MAX`.

    With ``detail=True`` returns ``{"floor", "contexts", "n_genes", "pooled", "clamped"}``.
    """
    if not 0.0 < q < 1.0:
        raise ValueError(f"context_floor: q must be in (0, 1), got {q!r}")
    keys: tuple = ()
    if tumor:
        codes, tis = resolve_context(tumor, path)
        keys = tuple(f"toil_tcga|{c}" for c in codes) or tuple(f"toil_gtex|{t}" for t in tis)
    elif tissue:
        _, tis = resolve_context(tissue, path)
        keys = tuple(f"toil_gtex|{t}" for t in tis)

    v, n = _floor_from(keys, q, path)
    pooled = not keys or v != v
    if v != v:
        _, ci, _, _ = _matrix(path)
        v, n = _floor_from(tuple(c for c in ci if c.startswith("toil_tcga|")), q, path)

    raw = v
    if prefilter and prefilter > 0:
        v = max(v, float(prefilter))
    if clamp:
        v = min(max(v, C_MIN), C_MAX)
    if not detail:
        return v
    return {"floor": v, "contexts": keys, "n_genes": n, "pooled": pooled, "clamped": v != raw}


def gene_level(gene: str, tumor: str | None = None, tissue: str | None = None,
               path: str | None = None) -> dict:
    """``{"tumor", "normal", "pan", "found"}`` -- one gene's level in TPM, three ways.

    All three from the single-pipeline table, so they are directly comparable and their differences
    mean something::

        tumor   the gene's median in that TCGA study, pooled over studies where the origin is
                an organ. ``None`` when no tumour type is given
        normal  the gene's median in the matched normal tissue(s). ``None`` when neither a tissue
                nor a tumour with a matched normal is given
        pan     the gene's median across the normal tissues where it is transcribed at all, and
                0.0 for a gene silent everywhere. **Always defined**, so a resolution chain has a
                last step that cannot fail

    ``found`` is whether the gene is in the reference at all. A gene that is not is not a gene with
    a level of zero, and the two must not be collapsed -- one is silence and the other is ignorance.
    """
    import numpy as np

    gi, ci, V, _ = _matrix(path)
    i = gi.get(str(gene).strip())
    if i is None:
        return {"tumor": None, "normal": None, "pan": None, "found": False}

    row = V[i]
    out: dict = {"found": True}

    tum_keys: tuple = ()
    nor_keys: tuple = ()
    if tumor:
        codes, tis = resolve_context(tumor, path)
        tum_keys = tuple(f"toil_tcga|{c}" for c in codes)
        nor_keys = tuple(f"toil_gtex|{t}" for t in tis)
    if tissue:
        _, tis = resolve_context(tissue, path)
        nor_keys = tuple(f"toil_gtex|{t}" for t in tis)

    def med(keys):
        cols = [ci[k] for k in keys if k in ci]
        return float(np.median(row[cols])) if cols else None

    out["tumor"] = med(tum_keys)
    out["normal"] = med(nor_keys)

    gcols = [j for c, j in ci.items() if c.startswith("toil_gtex|")]
    nz = row[gcols]
    nz = nz[nz > 0]
    out["pan"] = float(np.median(nz)) if nz.size else 0.0
    return out


def _f(x) -> float:
    """``float(x)`` or ``nan`` -- a submitted column may carry an empty string or ``None``."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def batch_scale(values, genes, tumor: str | None = None, path: str | None = None,
                detail: bool = False):
    """``(scale, n shared genes, fell back)`` -- what to divide a submitted column by to reach TPM.

    A median of ratios against the reference::

        scale = median over g of  x(g) / r(g)

    over the genes carrying a positive value in both, where ``r`` is the tumour type's own reference
    level, or the pooled reference where no tumour type is given. It is the standard robust
    per-sample size factor, and it does what a unit declaration cannot: it absorbs FPKM against TPM,
    raw counts against TPM, **and** one pipeline's TPM against another's, without the caller having
    to know which of those they have.

    Only genes positive in both count. A candidate at zero is a measurement of silence, not of
    scale, and letting it in drags the median toward zero in proportion to how many genes the tumour
    happens to have switched off.

    **Pass a whole-transcriptome profile, not a candidate list**, and the guards enforce it: the
    shared genes must clear :data:`MIN_SHARED` *and* cover :data:`MIN_COVERAGE` of the genes the
    reference context has switched on. A candidate list cannot, by two-and-a-half fold, and it must
    not -- a mutation is only called where it was seen in RNA, so candidate abundances are
    conditioned on detection and their ratio to the reference measures that conditioning. On the
    three screens carrying their own RNA-seq it returns 1.78, 2.18 and 3.15, all above 1, none of
    them a unit.

    Magnitude is deliberately not policed beyond :data:`GAMMA_MIN`--:data:`GAMMA_MAX`, which only
    catches a non-finite or structurally broken result: raw counts sit three or more orders of
    magnitude from TPM, and refusing that would reject the commonest case this exists to handle.

    ``detail=True`` adds the spread of the underlying ratios, ``(q75 - q25) / median``, and the
    coverage the estimate rests on. A pure rescale has spread 0; the measured screens run 2.1 to
    3.0. Tumour biology moves individual genes on its own, so spread is reported for a caller to
    judge rather than gated on a threshold -- coverage is the gate.

    >>> # a batch that is the reference times 7 recovers 7
    >>> batch_scale([7 * v for v in ref], names, tumor="SKCM")     # doctest: +SKIP
    (7.0, 412, False)
    """
    import numpy as np

    gi, ci, V, _ = _matrix(path)
    if tumor:
        codes, _t = resolve_context(tumor, path)
        keys = tuple(f"toil_tcga|{c}" for c in codes)
    else:
        keys = tuple(c for c in ci if c.startswith("toil_tcga|"))
    cols = [ci[k] for k in keys if k in ci]
    if not cols:
        return 1.0, 0, True

    # Vectorised deliberately. A whole-transcriptome profile is 20,000-30,000 genes, and the
    # readable version -- a Python loop with `np.median(V[i, cols])` inside it -- is one numpy call
    # per gene on a handful of values, which is all dispatch and no arithmetic. One gather and one
    # median along an axis does the same work in a single pass.
    idx = np.fromiter((gi.get(str(g).strip(), -1) for g in genes), dtype=np.int64,
                      count=len(genes))
    xv = np.asarray([_f(x) for x in values], dtype=float)
    ok = (idx >= 0) & np.isfinite(xv) & (xv > 0)
    if not ok.any():
        return (1.0, 0, True, float("nan"), 0.0) if detail else (1.0, 0, True)
    ref = np.median(V[np.ix_(idx[ok], cols)], axis=1)
    good = ref > 0
    num = xv[ok][good]
    den = ref[good]
    n = int(num.size)
    on = int((V[:, cols] > 0).any(axis=1).sum())
    cover = n / on if on else 0.0
    if n < MIN_SHARED or cover < MIN_COVERAGE:
        return (1.0, n, True, float("nan"), cover) if detail else (1.0, n, True)
    r = num / den
    s = float(np.median(r))
    spread = float((np.quantile(r, 0.75) - np.quantile(r, 0.25)) / s) if s else float("nan")
    if not np.isfinite(s) or not (GAMMA_MIN <= s <= GAMMA_MAX):
        return (1.0, n, True, spread, cover) if detail else (1.0, n, True)
    return (s, n, False, spread, cover) if detail else (s, n, False)
