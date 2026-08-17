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
import gzip
import os

__all__ = ["REFERENCE_FILE", "TUMOR_TISSUE", "matched_tissues",
           "fetch_reference", "load", "lookup", "impute", "tissues",
           "tumor_types", "safety_profile"]

#: Path inside the HF dataset repo.
REFERENCE_FILE = "expression/reference_expression.tsv.gz"

#: Columns of the reference table, in file order.
COLUMNS = ("key", "key_type", "source", "context", "median_tpm", "q25_tpm", "q75_tpm", "n")


def fetch_reference(path: str | None = None) -> str:
    """Local path to the reference table, downloading it from HF on first use.

    ``$MHCMATCH_EXPRESSION`` overrides, for offline and cluster runs, matching how
    ``$MHCMATCH_PMHC`` overrides the presentation panel."""
    if path:
        return path
    env = os.environ.get("MHCMATCH_EXPRESSION")
    if env:
        return env if os.path.isfile(env) else os.path.join(env, os.path.basename(REFERENCE_FILE))
    from huggingface_hub import hf_hub_download

    from .store import PMHC_REPO
    return hf_hub_download(repo_id=PMHC_REPO, repo_type="dataset", filename=REFERENCE_FILE)


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


def safety_profile(gene: str, top: int = 10, path: str | None = None) -> list[tuple[str, float]]:
    """``[(tissue, median_tpm)]`` for a gene across normal tissues, highest first.

    The safety read: a target expressed only in the tumour's lineage is very different from one
    expressed in heart and lung too, and the ranking score alone does not show that."""
    return _by_gene(fetch_reference(path)).get(gene, [])[:top]
