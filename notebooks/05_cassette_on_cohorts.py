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
    # 5 — Choosing what goes in the cassette, on a real trial's units

    **What this demonstrates.** Ranking gives you an ordered list. A construct needs a **set**, and
    the best set is not the top of the list. This notebook runs the whole applied chain from the shell
    on units that were actually manufactured and injected:

    ```
    rank pairs  ->  cassette select  ->  cassette score
    ```

    **What you should conclude.** Two things a ranking cannot tell you:

    * **Composition is not ranking.** Two cassettes with the same expected number of responding units
      can differ substantially in *P(at least one works)*, because the units in one of them fail
      together. The objective prices that; a sort cannot see it.
    * **`cassette score` is a cohort step, and that is not a detail.** Fit the calibration per donor
      and every donor's mean candidate probability becomes the declared prevalence by construction —
      so two donors' numbers are the same number and any cross-donor comparison built on them reads
      noise. Demonstrated below, not asserted.

    Data: `vaccines/ivac_mutanome_units.parquet` from `isalgo/pmhc_data` — the 125 units of the
    IVAC MUTANOME trial across 13 patients, with the published per-unit response calls. ~16 kB.

    **One stated assumption.** That deposit publishes no HLA typing, so the units are scored against a
    fixed common class-I panel named in the code below. It is an assumption, it is visible, and it is
    the reason this notebook does not report accuracy — it reports set geometry, which is what it is
    about.
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

    work = Path(tempfile.mkdtemp(prefix="mhcmatch-nb5-"))

    # **Stated, not typed from the trial**: IVAC MUTANOME publishes no HLA. A frequent European
    # class-I panel stands in, so the allotype channel has something real to work with. Every
    # allotype-flavoured number below is conditional on this choice and would move with a donor's
    # actual genotype -- which is exactly why `cassette select --universe` exists.
    PANEL = ["A*02:01", "A*01:01", "B*07:02", "B*08:01", "C*07:01", "C*07:02"]

    units = pl.read_parquet(fetch_file("vaccines/ivac_mutanome_units.parquet"))
    print(f"{units.height} manufactured units, {units['patient'].n_unique()} patients")
    print(f"units with a published response: {units['responded'].sum()} "
          f"({units['responded'].mean():.0%})")
    print(f"unit length: {units['peptide_len'].min()}-{units['peptide_len'].max()} aa "
          "(a vaccine unit is the long window, never the minimal epitope)")
    return PANEL, Path, np, pl, subprocess, sys, time, units, work


@app.cell
def _(mo):
    mo.md(r"""
    ## 5.1 A unit is a 27-mer; the model scores 9-mers

    A vaccine unit is the ~27-residue window around the mutation. Injecting the minimal epitope
    instead is not a smaller version of the right thing — a 9-mer loads onto any cell without
    costimulation and is the **tolerising** configuration — so the construct carries the window and
    the *scoring* happens on the epitopes inside it.

    Every 9-mer of every unit, against every allotype in the panel, then the best-scoring pair stands
    for the unit.
    """)
    return


@app.cell
def _(PANEL, pl, units, work):
    rows = []
    for r in units.iter_rows(named=True):
        p = r["peptide"]
        for i in range(len(p) - 8):
            for a in PANEL:
                rows.append({"unit_id": f'{r["patient"]}:{r["unit_id"]}', "patient": r["patient"],
                             "peptide": p[i:i + 9], "allele": a, "gene": r["gene"],
                             "tpm": r["exon_expression"], "responded": r["responded"],
                             "resp_cd8": r["resp_cd8"]})
    pairs = pl.DataFrame(rows)
    src = work / "pairs.tsv"
    pairs.write_csv(src, separator="\t")
    print(f"{pairs.height:,} (peptide, allele) pairs "
          f"= {units.height} units x ~{len(units['peptide'][0]) - 8} 9-mers x {len(PANEL)} allotypes")
    return pairs, src


@app.cell
def _(src, subprocess, sys, time, work):
    ranked_path = work / "ranked.tsv"
    _cmd = [sys.executable, "-m", "mhcmatch.cli", "rank", "pairs", str(src), "--cls", "mhc1",
            "--passthrough", "--prefix", "mm_", "--out", str(ranked_path)]
    print("$ mhcmatch " + " ".join(_cmd[3:]))
    _t0 = time.time()
    _p = subprocess.run(_cmd, capture_output=True, text=True)
    assert _p.returncode == 0, _p.stderr
    print(f"\nscored in {time.time() - _t0:.1f} s")
    print(_p.stderr.rstrip())
    return (ranked_path,)


@app.cell
def _(pl, ranked_path, work):
    _r = pl.read_csv(ranked_path, separator="\t", infer_schema_length=5000)
    # Select BEFORE renaming: the passthrough carries the caller's own `allele` (the panel member the
    # pair was built with) alongside our `mm_allele_scored` (the one that actually presented best),
    # and renaming into an occupied name is a duplicate-column error. Ours is the one that matters.
    pool = (_r.sort("mm_score", descending=True)
              .unique(subset=["unit_id"], keep="first")
              .select("patient", "unit_id", "peptide", "mm_allele_scored", "gene", "mm_score",
                      "responded", "resp_cd8")
              .rename({"patient": "donor", "mm_allele_scored": "allele", "mm_score": "score"})
              .sort("donor", "unit_id"))
    pool_path = work / "pool.tsv"
    pool.write_csv(pool_path, separator="\t")
    print(f"{pool.height} units, one row each, best epitope and its allotype")
    print(f"{pool['donor'].n_unique()} donors, {pool['allele'].n_unique()} allotypes used")
    print()
    print(pool.head(4))
    return pool, pool_path


@app.cell
def _(mo):
    mo.md(r"""
    ## 5.2 `cassette select` against ranking the same pool

    `select` maximises a certainty-equivalent objective: expected responding units, minus a risk term
    that charges for units failing *together*. Units share a failure when they share an allotype, when
    they look alike as sequence, or when one dominates another on score.

    A sort maximises the first term and is blind to the second. Below, both rules on the same pool,
    same *k*, same offset — and then the set properties each one bought.
    """)
    return


@app.cell
def _(pool_path, subprocess, sys, work):
    K = 8
    units_path = work / "units.tsv"
    _p = subprocess.run(
        [sys.executable, "-m", "mhcmatch.cli", "cassette", "select",
         "--candidates", str(pool_path), "-k", str(K), "--passthrough",
         "--out", str(units_path)],
        capture_output=True, text=True)
    assert _p.returncode == 0, _p.stderr
    print(_p.stderr.rstrip())
    return K, units_path


@app.cell
def _(K, np, pl, pool, units_path):
    chosen = pl.read_csv(units_path, separator="\t", infer_schema_length=2000)

    # The counterfactual: the same k, taken off the top of each donor's own ranking.
    top = (pool.sort("score", descending=True).group_by("donor").head(K))

    print(f"{'donor':>6} {'k':>3} {'designed':>28} {'top-k by score':>28}")
    print(f"{'':>6} {'':>3} {'allotypes  genes  caught':>28} {'allotypes  genes  caught':>28}")
    _rows = []
    for _d in sorted(pool["donor"].unique().to_list()):
        _c = chosen.filter(pl.col("donor") == _d)
        _t = top.filter(pl.col("donor") == _d)
        _stat = lambda f, col: f[col].n_unique()
        _row = (_d, _c.height,
                _stat(_c, "allele"), _stat(_c, "gene_in"), int(_c["responded"].sum()),
                _stat(_t, "allele"), _stat(_t, "gene"), int(_t["responded"].sum()))
        _rows.append(_row)
        print(f"{_row[0]:>6} {_row[1]:3d} {_row[2]:11d} {_row[3]:6d} {_row[4]:7d} "
              f"{_row[5]:15d} {_row[6]:6d} {_row[7]:7d}")

    _a = np.array([[r[2], r[3], r[4]] for r in _rows])
    _b = np.array([[r[5], r[6], r[7]] for r in _rows])
    print()
    print(f"{'totals':>6} {'':>3} {_a[:, 0].sum():11d} {_a[:, 1].sum():6d} {_a[:, 2].sum():7d} "
          f"{_b[:, 0].sum():15d} {_b[:, 1].sum():6d} {_b[:, 2].sum():7d}")
    print()
    print("`allotypes` and `genes` are what the design is spending capacity on. `caught` is how many")
    print("of the units the trial actually saw a response to are in each set -- reported because it")
    print("is the honest read-out, not because the design is optimised for it.")
    return chosen, top


@app.cell
def _(mo):
    mo.md(r"""
    ## 5.3 Why `cassette score` collects every donor first

    `rank` anchors `p_response` on the batch it is handed. Score each donor's cassette on its own and
    the offset is fitted to that donor's pool — so **every donor's mean candidate probability becomes
    the declared prevalence**, whatever their pool actually holds. Two donors' numbers are then the
    same number.

    Run it both ways on the same 13 cassettes:
    """)
    return


@app.cell
def _(np, pl, pool_path, subprocess, sys, units_path, work):
    def score(per_donor):
        o = work / f"scored_{'per' if per_donor else 'shared'}.tsv"
        cmd = [sys.executable, "-m", "mhcmatch.cli", "cassette", "score",
               "--cassettes", str(units_path), "--pool", str(pool_path), "--out", str(o)]
        if per_donor:
            cmd.append("--per-donor-offset")
        p = subprocess.run(cmd, capture_output=True, text=True)
        assert p.returncode == 0, p.stderr
        return pl.read_csv(o, separator="\t", infer_schema_length=2000)

    shared, per = score(False), score(True)
    print(f"{'':>22} {'offsets':>10} {'mean p_mean':>13} {'sd of p_mean':>14} {'spread of yield':>17}")
    for _name, _d in (("one over the cohort", shared), ("one per donor", per)):
        _off, _pm, _y = (_d["offset"].to_numpy(), _d["p_mean"].to_numpy(), _d["yield"].to_numpy())
        print(f"{_name:>22} {len(np.unique(np.round(_off, 9))):10d} {_pm.mean():13.6f} "
              f"{_pm.std():14.3e} {_y.max() - _y.min():17.4f}")
    _sh = shared["p_mean"].to_numpy().std()
    _pd = per["p_mean"].to_numpy().std()
    _ys = shared["yield"].to_numpy()
    _yp = per["yield"].to_numpy()
    print()
    print(f"per-donor calibration shrinks the spread of `p_mean` across donors "
          f"{_sh / max(_pd, 1e-30):.1f}x")
    print(f"and the spread of `yield` {(_ys.max() - _ys.min()) / max(_yp.max() - _yp.min(), 1e-30):.1f}x.")
    print()
    print("The direction is the point and the size depends on the pool. These donors have ~10")
    print("candidates each, so `prob_offset` cannot land every mean exactly on the prevalence. On a")
    print("real cohort it does: the benchmark's 7,261 TCGA donors, each with hundreds to thousands of")
    print("candidates, all land on a pool mean of 0.060163 with a standard deviation of 2.75e-17 --")
    print("thirteen orders of magnitude below the spread the cohort offset keeps. At that point two")
    print("donors' numbers are literally the same number and a triage built on them reads noise.")
    return per, score, shared


@app.cell
def _(mo):
    mo.md(r"""
    ## 5.4 `lam` is relative by construction, and that is a different guarantee

    `yield` is a **level**: how many units this cassette is expected to get a response from. A level
    only means something against a fixed calibration, which is why the cohort offset above matters.

    `lam` is a **contrast**: nats above a uniform random subset of *that donor's own pool*. Both terms
    move together when the offset moves, so `lam` stays a statement about this donor's selection
    against this donor's background under either calibration — which is what makes it quotable when
    two cassettes were not drawn from one pool.

    It is **not** invariant to the offset, and on pools this small it visibly is not. Below: the two
    calibrations side by side, and what actually survives.
    """)
    return


@app.cell
def _(np, per, shared):
    _a = shared.sort("donor")
    _b = per.sort("donor")
    print(f"{'donor':>6} {'k':>3} {'lam (cohort offset)':>21} {'lam (per donor)':>18}")
    for _r1, _r2 in zip(_a.iter_rows(named=True), _b.iter_rows(named=True)):
        print(f"{_r1['donor']:>6} {_r1['k']:3d} {_r1['lam']:21.4f} {_r2['lam']:18.4f}")
    _x, _y = _a["lam"].to_numpy(), _b["lam"].to_numpy()
    print()
    print(f"largest absolute difference between the two calibrations: {np.abs(_x - _y).max():.4f}")
    print(f"Spearman rho between them across the {len(_x)} donors: "
          f"{np.corrcoef(np.argsort(np.argsort(_x)), np.argsort(np.argsort(_y)))[0, 1]:+.3f}")
    print()
    print("So `lam` MOVES with the calibration -- it is not an invariant, and a notebook that claimed")
    print("it was would be wrong. What survives is that it stays positive and stays a within-donor")
    print("contrast: under both calibrations every cassette reads above a random subset of its own")
    print("donor's pool, which is the claim `lam` actually makes. `yield`, by contrast, stops being")
    print("comparable between donors at all once the offset is fitted per donor.")
    print()
    print(f"cassettes reading above their own donor's random subset: "
          f"{int((_x > 0).sum())}/{len(_x)} (cohort offset), "
          f"{int((_y > 0).sum())}/{len(_y)} (per donor)")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to take away

    * A cassette is a **set**, and the top of a ranking is not the best set. What a sort cannot see is
      units failing together — same allotype, same gene, one dominating another.
    * **`cassette score` fits one offset over the whole run.** Per donor it reports an *enrichment*
      against that donor's own background, which is a real quantity and no longer comparable between
      donors. The flag exists so the choice is deliberate; the default is the comparable one.
    * **`lam` is a contrast, `yield` is a level.** `lam` still moves with the offset — it is not an
      invariant — but it stays a statement about this donor's selection against this donor's own
      background, which is what makes it the one to quote when two cassettes came from different
      pools. `yield` stops being comparable between donors the moment the offset is fitted per donor.
    * A vaccine unit is the **long window**. The minimal epitope is the tolerising configuration, and
      the library will not silently build one from a `peptide` column.

    Next: [notebook 6](06_assembly_and_safety.py) — turning the chosen units into a construct.
    """)
    return


if __name__ == "__main__":
    app.run()
