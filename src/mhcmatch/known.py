"""Known-epitope reference sets, assembled from the public deposits, for exact-match lookup.

**An exact match is stronger evidence than any model output, so it is a flag and never a score.**
:mod:`mhcmatch.rank` reports which set a candidate was found in and floats those candidates into a
tier of their own, with the model score still shown beside them. Burying "this peptide is a
confirmed NCI neoantigen" inside a weighted sum lets a mediocre model score dilute the one piece of
direct evidence in the row.

Five sets, each answering a different question about a candidate:

===================  ======================================================================
set                  a hit means
===================  ======================================================================
``neoantigen``       a **confirmed immunogenic tumour neoantigen** -- the NCI/Gartner
                     deconvolved minimal peptides, the epitope-resolution screens (NCI,
                     HiTIDE, TESLA) and the aggregated cohorts, all restricted to rows the
                     assay called positive. This candidate has already been shown to work.
``neoantigen_neg``   **screened and found non-immunogenic.** The same deposits, negative
                     rows. Not evidence of nothing: it is the one label that says this exact
                     peptide was tested and did not respond.
``immunogenic``      a positive T-cell assay in IEDB against any source -- viral, bacterial,
                     tumour. Immunogenic in some context, not necessarily this one.
``self``             present in the **thymic self-immunopeptidome** (HLA Ligand Atlas). The
                     tolerance argument, and a cross-reactivity/autoimmunity flag for a
                     vaccine: reactive T cells were plausibly deleted.
``viral``            a presented pathogen-derived peptide. A pre-existing anti-pathogen
                     repertoire may cross-react.
===================  ======================================================================

``neoantigen`` and ``neoantigen_neg`` are disjoint by construction: a peptide reported positive
anywhere is a positive, because a negative call in one assay does not overturn a positive one in
another. That rule is applied **after** pooling, in :func:`load`.

Every set is fetched from the public HF dataset (:data:`mhcmatch.store.PMHC_REPO`) and cached, so a
fresh install needs no pre-staged data. ``MHCMATCH_PMHC_DIR`` points at a local mirror instead.

    >>> from mhcmatch import known
    >>> refs = known.load()                                  # doctest: +SKIP
    >>> "GILGFVFTL" in refs["viral"]                         # doctest: +SKIP
    True
"""
from __future__ import annotations

import csv
import functools
import gzip

__all__ = ["SOURCES", "SET_NAMES", "load", "lookup"]

#: ``set -> [(repo-relative file, label column or None, values that count as a hit)]``. A ``None``
#: label column means every row of the file belongs to the set.
SOURCES: dict[str, list[tuple]] = {
    "neoantigen": [
        ("neoantigens/nci_gartner_mmp.tsv.gz", "immunogenicity", {"cd8", "1", "positive"}),
        ("neoantigens/neoantigens_tested_peptides.tsv.gz", "immunogenicity",
         {"cd8", "1", "positive"}),
        ("neoantigens/neoag_tested.tsv.gz", "immunogenicity", {"1", "cd8", "positive"}),
        ("neoantigens/neoag_tested_hsa.tsv.gz", "immunogenicity", {"1", "cd8", "positive"}),
        ("neoantigens/neoag_tested_mmu.tsv.gz", "immunogenicity", {"1", "cd8", "positive"}),
    ],
    "neoantigen_neg": [
        ("neoantigens/neoantigens_tested_peptides.tsv.gz", "immunogenicity",
         {"0", "negative", "non-immunogenic"}),
        ("neoantigens/neoag_tested.tsv.gz", "immunogenicity", {"0", "negative"}),
        ("neoantigens/neoag_tested_hsa.tsv.gz", "immunogenicity", {"0", "negative"}),
        ("neoantigens/neoag_tested_mmu.tsv.gz", "immunogenicity", {"0", "negative"}),
    ],
    "immunogenic": [("immunogenicity/chowell_rebuilt.tsv.gz", "label", {"1"})],
    "self": [("thymus/thymus_immunopeptidome.tsv.gz", None, None)],
    "viral": [("ligandome/viral_foreign_iedb.tsv.gz", None, None)],
}
#: Report order: strongest direct evidence first, so :func:`lookup` names the best one.
SET_NAMES = ("neoantigen", "neoantigen_neg", "immunogenic", "self", "viral")


def _read(rel: str, label_col: str | None, hit_values: set | None) -> set[str]:
    from .store import fetch_file

    path = fetch_file(rel)
    op = gzip.open if path.endswith(".gz") else open
    out: set[str] = set()
    with op(path, "rt", encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            p = (row.get("peptide") or "").strip().upper()
            if not p or not p.isalpha():
                continue
            if label_col is not None:
                v = (row.get(label_col) or "").strip().lower()
                if v not in hit_values:
                    continue
            out.add(p)
    return out


@functools.lru_cache(maxsize=4)
def load(names: tuple[str, ...] | None = None) -> dict[str, frozenset]:
    """``{set name: frozenset of peptides}``, downloaded on first use and cached.

    ``names`` restricts which sets are built -- each one costs a download and a full-file scan, so
    ask for what you need. Default is all of :data:`SET_NAMES`.

    A peptide reported immunogenic anywhere is removed from ``neoantigen_neg``: one assay calling it
    negative does not overturn another calling it positive, and leaving it in both would make the
    flag report whichever set happened to be checked first."""
    want = tuple(names) if names else SET_NAMES
    bad = [n for n in want if n not in SOURCES]
    if bad:
        raise ValueError(f"unknown reference set(s) {bad}; expected from {SET_NAMES}")
    out = {n: set().union(*(_read(*spec) for spec in SOURCES[n])) for n in want}
    if "neoantigen_neg" in out and "neoantigen" in out:
        out["neoantigen_neg"] -= out["neoantigen"]
    return {n: frozenset(v) for n, v in out.items()}


def lookup(peptide: str, refs: dict | None = None) -> str:
    """The **first** set in :data:`SET_NAMES` containing this peptide exactly, else ``""``.

    Order is deliberate: a peptide that is both a confirmed neoantigen and a thymic self-peptide
    should report as the neoantigen, because that is the stronger and more actionable evidence --
    the self hit is still visible by checking ``refs`` directly."""
    p = peptide.strip().upper()
    r = refs if refs is not None else load()
    for name in SET_NAMES:
        if p in r.get(name, ()):
            return name
    return ""


def demo() -> None:
    """Self-check: run with ``python -m mhcmatch.known`` (downloads ~5 MB on first use)."""
    refs = load()
    for n in SET_NAMES:
        assert n in refs and len(refs[n]) > 0, n
    # the two neoantigen sets are disjoint after the positive-wins rule
    assert not (refs["neoantigen"] & refs["neoantigen_neg"])
    # every peptide is a canonical sequence
    assert all(p.isalpha() and p.isupper() for p in list(refs["neoantigen"])[:200])
    # lookup respects the report order, and an unknown peptide reports nothing
    a_neo = next(iter(refs["neoantigen"]))
    assert lookup(a_neo, refs) == "neoantigen"
    assert lookup(a_neo.lower(), refs) == "neoantigen"          # case-insensitive
    assert lookup("WWWWWWWWWWWW", refs) == ""
    # a set can be built on its own without paying for the rest
    assert set(load(("viral",))) == {"viral"}
    print("ok - " + ", ".join(f"{n} {len(refs[n]):,}" for n in SET_NAMES))


if __name__ == "__main__":
    demo()
