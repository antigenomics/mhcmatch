import marimo

__generated_with = "0.23.16"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md(r"""
    # 4 — Molecular mimicry and self (`mhcmatch.mimics`)

    **What this demonstrates.** For a strong binder, `mimics.scan` searches reference peptide sets for
    **mimics** — near-identical presented peptides — and reports a presentation-aware E-value per
    category. Three references ship as `DEFAULT_REFS`:

    | category | reference | what a hit means |
    |---|---|---|
    | `thymus` | thymic self-immunopeptidome (HLA Ligand Atlas) | the peptide resembles something presented during **negative selection** — reactive T cells were likely deleted, *and* it flags cross-reactivity / autoimmune risk for a vaccine |
    | `viral` | foreign presented peptides | a pre-existing anti-pathogen repertoire may cross-react — mimicry that can *raise* immunogenicity |
    | `neoag` | the tested-neoantigen database | has this (or something near-identical) been reported before |

    **What you should conclude.** Read section 4.3 before you build any pathogen-similarity feature.
    Three canonical epitopes — `GILGFVFTL`, `NLVPMVATV` and `KLGGALQAK` — all return `n_exact = 1`
    against the viral reference **because they are in it**. That is circularity reported as signal. A
    mimicry feature must exclude the query's own identity (and ideally its source study) from the
    reference, or it measures nothing but its own input.

    `find_mimics` already excludes the exact query inside the fuzzy search; the `n_exact`
    set-membership check in `scan` does not, and that asymmetry is what section 4.3 makes visible.

    All three references bootstrap from the public HuggingFace dataset `isalgo/pmhc_data`.
    """)
    return


@app.cell
def _():
    import os
    import time

    from huggingface_hub import hf_hub_download

    from mhcmatch import mimics
    from mhcmatch.store import PMHC_REPO

    CLS, SPECIES = "mhc1", "human"

    _t0 = time.time()
    # Fetch each DEFAULT_REFS file from the public HF dataset. They land in one cached snapshot
    # directory, whose root is exactly the `pmhc_dir` that load_reference_sets expects.
    _paths = [
        hf_hub_download(repo_id=PMHC_REPO, repo_type="dataset", filename=_rel)
        for _rel, _kind in mimics.DEFAULT_REFS.values()
    ]
    pmhc_dir = os.path.commonpath(_paths)
    print(f"references fetched/cached in {time.time() - _t0:.1f} s")
    for _rel, _kind in mimics.DEFAULT_REFS.items():
        print(f"  {_rel:<8s} {_kind[0]:<44s} kind={_kind[1]}")

    _t0 = time.time()
    self_set, foreign_sets = mimics.load_reference_sets(pmhc_dir, CLS, SPECIES)
    print(f"\nload_reference_sets: {time.time() - _t0:.1f} s")
    print(f"  thymus (self) {len(self_set):>7d} peptides")
    for _k, _v in foreign_sets.items():
        print(f"  {_k:<13s} {len(_v):>7d} peptides")
    return CLS, foreign_sets, mimics, self_set, time


@app.cell
def _(mo):
    mo.md(r"""
    ## 4.1 Scanning binders

    `scan` takes `(peptide, allele)` pairs and makes **one** `find_mimics` call per binder, which
    scores every category at once. It returns one `MimicResult` per (binder, category) that has at
    least one same-length reference peptide within `near_subs` substitutions — categories with nothing
    to report are simply absent from the output.
    """)
    return


@app.cell
def _(CLS, foreign_sets, mimics, mo, self_set, time):
    BINDERS = [
        ("GILGFVFTL", "HLA-A*02:01"),  # influenza A M1
        ("NLVPMVATV", "HLA-A*02:01"),  # CMV pp65
        ("KLGGALQAK", "HLA-A*03:01"),  # CMV IE1
        ("NLVPMVATL", "HLA-A*02:01"),  # a single-substitution variant -- NOT in any reference
    ]

    _t0 = time.time()
    results = mimics.scan(BINDERS, self_set, foreign_sets, cls=CLS)
    print(f"scan of {len(BINDERS)} binders over 3 references: {time.time() - _t0:.1f} s")
    print(f"{len(results)} (binder, category) results")

    def esc(allele):
        """Escape the `*` in an HLA name so markdown does not read it as emphasis."""
        return allele.replace("*", r"\*")


    _rows = "\n".join(
        f"| `{_r.binder}` | {esc(_r.allele)} | {_r.category} | {_r.n_exact} | "
        f"{_r.n_near} | `{_r.top_mimic}` | {_r.top_subs} | {_r.e_value:.3g} | {_r.n_hits} |"
        for _r in results
    )
    mo.md(
        f"""
    | binder | allele | category | `n_exact` | `n_near` | `top_mimic` | `top_subs` | `e_value` | `n_hits` |
    |---|---|---|--:|--:|---|--:|--:|--:|
    {_rows}
    """
    )
    return BINDERS, results


@app.cell
def _(BINDERS, mimics, mo, results):
    summary = mimics.patient_summary(results, BINDERS)
    _rows = "\n".join(f"| `{_k}` | {_v} |" for _k, _v in summary.items())
    mo.md(
        f"""
    ## 4.2 Patient-level roll-up

    `patient_summary` aggregates into the counts a dashboard row needs. It takes the **full** binder
    list as well, so binders with zero mimics are still counted in the denominator.

    | field | value |
    |---|--:|
    {_rows}

    `n_tolerance_risk` counts binders with any significant thymic/self mimic — the cross-reactivity and
    autoimmunity flag. `mimics.write_table(results, path)` writes the per-(binder, category) rows as a
    TSV with the `NATIVE_COLUMNS` schema.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 4.3 The leakage trap

    Look again at the `n_exact` column in section 4.1.

    `GILGFVFTL`, `NLVPMVATV` and `KLGGALQAK` each score `n_exact = 1` against the viral reference — and
    against the tested-neoantigen reference too. Not because they resemble a pathogen peptide, but
    because **they are canonical, heavily-published epitopes and the reference sets are built from
    published data**. Each of them is matching itself.

    The fourth query, `NLVPMVATL`, is not in any reference, and it behaves the way a real query should:
    `n_exact = 0`, with genuine near mimics one substitution away.

    Re-run the same scan with the queries removed from every reference and the difference is stark.
    """)
    return


@app.cell
def _(BINDERS, CLS, foreign_sets, mimics, mo, results, self_set, time):
    _queries = {_p for _p, _a in BINDERS}
    clean_self = [_p for _p in self_set if _p not in _queries]
    clean_foreign = {_k: [_p for _p in _v if _p not in _queries] for _k, _v in foreign_sets.items()}
    print(f"removed {len(self_set) - len(clean_self)} query peptides from thymus")
    for _k in foreign_sets:
        print(f"removed {len(foreign_sets[_k]) - len(clean_foreign[_k])} query peptides from {_k}")

    _t0 = time.time()
    clean_results = mimics.scan(BINDERS, clean_self, clean_foreign, cls=CLS)
    print(f"\nde-leaked scan: {time.time() - _t0:.1f} s -> {len(clean_results)} results "
          f"(was {len(results)})")

    _before = {(_r.binder, _r.category): _r for _r in results}
    _after = {(_r.binder, _r.category): _r for _r in clean_results}
    _keys = sorted(set(_before) | set(_after))
    _rows = "\n".join(
        f"| `{_b}` | {_c} | "
        f"{_before[(_b, _c)].n_exact if (_b, _c) in _before else '-'} | "
        f"{_after[(_b, _c)].n_exact if (_b, _c) in _after else '-'} | "
        f"{_before[(_b, _c)].n_near if (_b, _c) in _before else '-'} | "
        f"{_after[(_b, _c)].n_near if (_b, _c) in _after else '-'} | "
        f"{('`' + _after[(_b, _c)].top_mimic + '`') if (_b, _c) in _after else 'dropped'} |"
        for _b, _c in _keys
    )
    mo.md(
        f"""
    | binder | category | `n_exact` before | `n_exact` after | `n_near` before | `n_near` after | `top_mimic` after |
    |---|---|--:|--:|--:|--:|---|
    {_rows}

    ### Read this row by row

    * Every `n_exact = 1` collapses to `0`. All of it was self-identity.
    * `NLVPMVATV` and `KLGGALQAK` **disappear from the viral category entirely** — they had no near
      mimics at all, so their whole "viral similarity" was the leak.
    * `GILGFVFTL` keeps 3 genuine near mimics; the top one is one substitution away, which is a real
      cross-reactivity signal.
    * `NLVPMVATL` is unchanged, because it was never in the reference to begin with.

    **The rule.** Any pathogen- or self-similarity feature must exclude the query's own identity from
    the reference before it is computed. Ideally exclude the query's source study too: a peptide
    published once appears in several derived compendia, so identity-only filtering can still leave a
    near-duplicate behind.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---

    ### Notes

    * **`scan` scores cross-reactivity, not presentation or immunogenicity.** Compose it with the
      presentation / affinity / binder scores from `mhcmatch.predict` (notebook 1) rather than reading
      it on its own.
    * **`refs=` overrides `DEFAULT_REFS`** if you want a different reference mix — the value is
      `{name: (path_under_pmhc_dir, kind)}`, where exactly one entry carries `kind="self"` and becomes
      the tolerance reference.
    * `max_subs` is the fuzzy-search radius; `near_subs` is the reporting radius for `n_near`. They are
      separate on purpose — search wide, report narrow.
    * The `e_value` column is the raw presentation-aware search statistic from
      `mhcmatch.search.find_mimics`, kept for reference; it scales with the size of the reference set,
      so compare it within a category, never across.
    """)
    return


if __name__ == "__main__":
    app.run()
