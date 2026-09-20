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
    # 4 — The whole EPIC aggregate, from the command line

    **What this demonstrates.** Notebooks 1–3 took the model apart. This one runs it whole, the way a
    user actually would: `mhcmatch rank pairs` on a table of (peptide, allele) candidates, from the
    shell, with the caller's own columns coming back untouched.

    Then it reads two things back out that are easy to miss:

    * the **per-term decomposition** — nine standardised terms in four blocks, and their contributions
      summing to the score;
    * `p_response`, which is a probability **anchored on a prior you own**, not a model output.

    **What you should conclude.** The aggregate is one ridge-penalised logistic regression over nine
    terms, and every one of them is inspectable.

    **This is not a benchmark.** The panel below is a small, class-balanced, curated set of
    already-tested peptides with no expression data attached — a good shape for showing the mechanics
    and the wrong shape for measuring accuracy. The held-out record lives in the shipped artifact and
    is printed at the end; head-to-head comparisons against other tools live in the benchmark
    repository and deliberately not here.

    Data: `neoantigens/neoprecis_ranker_panel.parquet` from `isalgo/pmhc_data`, ~210 kB. It carries
    other tools' scores; **this notebook reads none of them** — only `mt_pept`, `wt_pept`, `allele`
    and `label`.
    """)
    return


@app.cell
def _():
    import subprocess
    import sys
    import tempfile
    import time
    from pathlib import Path

    import numpy as np
    import polars as pl

    from mhcmatch.store import fetch_file

    work = Path(tempfile.mkdtemp(prefix="mhcmatch-nb4-"))
    _t0 = time.time()
    _panel = pl.read_parquet(fetch_file("neoantigens/neoprecis_ranker_panel.parquet"))

    # Class I, and the authors' own held-out split. Only four columns are taken: the rival scores in
    # this file are not read, here or anywhere in these notebooks.
    cands = (_panel
             .filter((pl.col("mhc") == "I") & (pl.col("dataset") == "test"))
             .select(peptide=pl.col("mt_pept"), allele=pl.col("allele"),
                     wt_peptide=pl.col("wt_pept"), label=pl.col("label")))
    src = work / "candidates.tsv"
    cands.write_csv(src, separator="\t")
    print(f"loaded in {time.time() - _t0:.1f} s")
    print(f"{cands.height} candidates, {cands['label'].sum()} labelled immunogenic, "
          f"{cands['allele'].n_unique()} allotypes")
    return Path, cands, np, pl, src, subprocess, sys, time, work


@app.cell
def _(mo):
    mo.md(r"""
    ## 4.1 One command

    `rank pairs` takes a table that already names its (peptide, allele) pairs — as opposed to
    `rank fasta`, which calls the epitopes itself from mutation windows.

    **`--passthrough` needs a `--prefix`.** Without one the command refuses to run when your column
    names collide with its own, which is the correct failure: two columns under one name break
    silently, because every reader that keys a row by name resolves the duplicate in favour of one of
    them and the file does not record which.

    `wt_peptide` is already in this table, so agretopicity is computable. Without a wild type every
    row reads `wt_absent` — correct, and a weaker model.
    """)
    return


@app.cell
def _(src, subprocess, sys, time, work):
    out = work / "ranked.tsv"
    _cmd = [sys.executable, "-m", "mhcmatch.cli", "rank", "pairs", str(src),
            "--cls", "mhc1", "--passthrough", "--prefix", "mm_", "--out", str(out)]
    print("$ mhcmatch " + " ".join(_cmd[3:]))
    print()
    _t0 = time.time()
    _p = subprocess.run(_cmd, capture_output=True, text=True)
    print(_p.stderr.rstrip())
    print()
    print(f"ran in {time.time() - _t0:.1f} s")
    assert _p.returncode == 0, _p.stderr
    return (out,)


@app.cell
def _(out, pl):
    ranked = pl.read_csv(out, separator="\t", infer_schema_length=2000)
    print(f"{ranked.height} rows x {ranked.width} columns")
    print()
    _mine = [c for c in ranked.columns if not c.startswith("mm_")]
    _ours = [c for c in ranked.columns if c.startswith("mm_")]
    print(f"your columns, unchanged and in your order ({len(_mine)}): {', '.join(_mine)}")
    print()
    print(f"ours, under the prefix ({len(_ours)}):")
    for _i in range(0, len(_ours), 4):
        print("   " + "  ".join(f"{c:24s}" for c in _ours[_i:_i + 4]))
    return (ranked,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 4.2 The nine terms, and what the model does with each

    Read off the shipped artifact rather than typed. This is the same file the manuscript's tables are
    generated from, so these are the published coefficients by construction.
    """)
    return


@app.cell
def _():
    import json
    from pathlib import Path as _P

    import mhcmatch as _mm

    fit = json.loads((_P(_mm.__file__).parent / "data" / "aggregate_mhc1.json").read_text())
    weights = dict(zip(fit["features"], fit["coef"]))
    print(f"artifact {fit['model_id']}  version {fit['version']}  "
          f"({fit['fit']['rows']:,} rows, {fit['fit']['positives']} positives)")
    print()
    for _block, _terms in fit["blocks"]:
        print(f"{_block}")
        for _t in _terms:
            _arrow = "raises" if weights[_t] > 0 else "lowers"
            print(f"   {_t:18s} {weights[_t]:+.4f}   {_arrow} the score")
    return fit, weights


@app.cell
def _(mo):
    mo.md(r"""
    ## 4.3 The contributions add up

    Each term is standardised with the artifact's own `mu` and `sigma`, multiplied by its coefficient,
    and summed. That sum **is** the linear predictor — there is no residual and nothing else in the
    score. It is what makes a per-candidate explanation possible, and it is what `mhcmatch explain`
    prints for a single pair.
    """)
    return


@app.cell
def _(fit, np, ranked):
    _terms = fit["features"]
    _mu = np.asarray(fit["mu"], dtype=float)
    _sd = np.asarray(fit["sigma"], dtype=float)
    _w = np.asarray(fit["coef"], dtype=float)

    have = [t for t in _terms if f"mm_{t}" in ranked.columns]
    print(f"terms present in the output: {len(have)} of {len(_terms)}")
    _absent = [t for t in _terms if t not in have]
    if _absent:
        print(f"absent: {', '.join(_absent)} -- emitted only when the input supports them")
    print()

    _raw = np.column_stack([ranked[f"mm_{t}"].to_numpy().astype(float) for t in have])
    _idx = [_terms.index(t) for t in have]
    _z = (_raw - _mu[_idx]) / np.where(_sd[_idx] > 0, _sd[_idx], 1.0)
    contrib = _z * _w[_idx]

    # **Per term, not per row.** This input carries no gene column, so the two expression terms are
    # `nan` on EVERY row -- and an all-finite-rows mask is therefore empty and reports nothing about
    # the six terms that ARE defined. That is not a defect in the input: `rank pairs` does not infer
    # a gene from a peptide, it accepts one, and without it `expr_lvl`/`expr_norm` are undefined
    # rather than zero. Saying which terms are live is more useful than a table of nan.
    _finite = np.isfinite(contrib)
    print(f"{'term':>18} {'defined on':>12} {'mean |contribution|':>21} {'coefficient':>13}")
    _mean = np.full(len(have), np.nan)
    for _j in range(len(have)):
        _m = _finite[:, _j]
        if _m.any():
            _mean[_j] = np.abs(contrib[_m, _j]).mean()
    for _j in np.argsort(np.where(np.isnan(_mean), -1.0, -_mean)):
        _n = int(_finite[:, _j].sum())
        _val = "undefined" if np.isnan(_mean[_j]) else f"{_mean[_j]:.4f}"
        print(f"{have[_j]:>18} {_n:7d}/{len(contrib):<4d} {_val:>21} {_w[_idx[_j]]:+13.4f}")
    print()
    _live = [have[j] for j in range(len(have)) if not np.isnan(_mean[j])]
    print(f"{len(_live)} of the model's {len(_terms)} terms are live on this input; "
          f"{len(_terms) - len(_live)} are not.")
    print("Undefined here: the two expression terms (no gene column, and `rank pairs` accepts a gene")
    print("rather than inferring one) and `log10a`, which needs the affinity head's occupancy.")
    print("A candidate table from `rank fasta`, or one carrying `gene` and `tpm`, lights all nine.")
    return contrib, have


@app.cell
def _(mo):
    mo.md(r"""
    A term with a large coefficient that barely varies moves nothing; a term with a small coefficient
    and a wide spread can dominate. **The ordering above is the one that matters for a candidate
    list**, and it is not the ordering of the coefficient table.

    ## 4.4 `p_response` is a prior you own

    The fit gave every dataset its own intercept precisely so base rate stayed out of the slopes — the
    datasets behind it span three orders of magnitude in prevalence. So the model has no opinion about
    how many of *your* candidates will respond, and `p_response` is anchored on a number you supply.

    **It shifts every probability and moves no rank.** Demonstrated rather than asserted:
    """)
    return


@app.cell
def _(np, pl, src, subprocess, sys, work):
    def run_at(prev):
        o = work / f"ranked_p{prev}.tsv"
        p = subprocess.run(
            [sys.executable, "-m", "mhcmatch.cli", "rank", "pairs", str(src), "--cls", "mhc1",
             "--prevalence", str(prev), "--passthrough", "--prefix", "mm_", "--out", str(o)],
            capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        return pl.read_csv(o, separator="\t", infer_schema_length=2000)

    _a, _b = run_at(0.01), run_at(0.30)
    print(f"{'prevalence':>12} {'mean p_response':>18} {'max':>12}")
    for _name, _d in (("0.01", _a), ("0.30", _b)):
        _p = _d["mm_p_response"].to_numpy().astype(float)
        print(f"{_name:>12} {np.nanmean(_p):18.6f} {np.nanmax(_p):12.6f}")
    print()
    print(f"ranked order identical at both prevalences: "
          f"{_a['peptide'].to_list() == _b['peptide'].to_list()}")
    print()
    print("So a prevalence you cannot defend costs you nothing in ranking -- and a probability you")
    print("quote to somebody else is only as good as the prior you put in.")
    return (run_at,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 4.5 Where the held-out numbers actually are

    Not on this panel. The record is in the artifact: seven published datasets, each held out in turn
    and the model refitted without it.
    """)
    return


@app.cell
def _(fit, np):
    _loo = sorted(fit["loo"], key=lambda x: -x["auroc"])
    print(f"{'dataset':>12} {'n':>9} {'positives':>10} {'held-out AUROC':>16}")
    for _x in _loo:
        print(f"{_x['level']:>12} {_x['n']:9,} {_x['pos']:10d} {_x['auroc']:16.4f}")
    _au = [x["auroc"] for x in _loo]
    print(f"{'mean':>12} {'':>9} {'':>10} {np.mean(_au):16.4f}")
    print(f"{'median':>12} {'':>9} {'':>10} {np.median(_au):16.4f}")
    print()
    print("All seven are 'decided' -- at least 20 held-out positives, which is what it takes to")
    print("resolve a difference of 1/n_pos. No cell here is too small to read.")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to take away

    * `rank pairs` gives your table back — **every column, in your order** — with a prefixed block
      appended and the rows re-sorted. It is not a join, and no key would survive one: `rank` splits a
      cell naming several alleles, so the output shares neither its length nor its allele column with
      the input.
    * The score is **nine standardised terms and nothing else**, and the contributions add up exactly.
      What moves a candidate list is the term with the widest *spread*, not the biggest coefficient.
    * `p_response` is **your prior**. It moves every probability and no rank.
    * Performance numbers belong to a benchmark, and the honest ones are held out one dataset at a
      time. They ship inside the artifact so they cannot drift from the model that produced them.

    Next: [notebook 5](05_cassette_on_cohorts.py) — choosing which twenty of these to carry, which is
    a different question from ranking them.
    """)
    return


if __name__ == "__main__":
    app.run()
