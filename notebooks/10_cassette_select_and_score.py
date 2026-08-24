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
    # 10 — Choosing the units, and scoring the cassette (`mhcmatch.cassette`)

    **What this demonstrates.** The two operations on a *whole published corpus*, not a toy:
    `cassette select` on the 46 patients of the NCI gastrointestinal screen held out of the EPIC
    fit, and `cassette score` on the 1,631 units two trials actually manufactured. Notebook 9 shows
    *why* a cassette is a set problem; this one is the operation.

    **What you should conclude.**

    1. `select` maximises **mean minus variance** of the responding-unit count. A sort maximises the
       mean. They are different objectives and each wins on its own — reporting only one of them
       decides the comparison before making it.
    2. `sum p` is a **level** and it is comparable only if every cassette was calibrated together.
       Fit the offset per donor and every donor's mean becomes the declared prevalence, exactly.
    3. `lam` needs no shared calibration at all: it divides by the donor's own pool, so it compares
       cassettes across donors **and** across sizes.
    4. Selection is worth a great deal on an unfiltered pool and very little on a filtered one, and
       the difference is a base rate.

    Everything is fetched from the public `isalgo/pmhc_data` deposit, so this runs on a bare
    `pip install mhcmatch`. Roughly 90 s, most of it the corpus channels.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. A donor's unfiltered candidate list

    `nci_parkhurst_gi.parquet` is the NCI Surgery Branch gastrointestinal screen: every non-silent
    mutation the patient carried, tiled and assayed against autologous TIL. 29 of its 75 patients are
    inside the corpus the shipped EPIC scorer was fitted on, and the deposit ships `in_epic_fit` as a
    column so that cannot be missed. **We keep the 46 that are not.**
    """)
    return


@app.cell
def _():
    import numpy as np
    import polars as pl

    from mhcmatch import cassette as CA
    from mhcmatch import rank as RK
    from mhcmatch.store import fetch_file

    d = pl.read_parquet(fetch_file("neoantigens/nci_parkhurst_gi.parquet"))
    held = d.filter(~pl.col("in_epic_fit"))
    print(f"{d.height:,} screened mutations, {d['patient'].n_unique()} patients")
    print(f"held out of the EPIC fit: {held['patient'].n_unique()} patients, "
          f"{held.height:,} mutations, {int(held['cd8'].sum())} CD8+, {int(held['cd4'].sum())} CD4+")
    pool_sizes = held.group_by("patient").len()["len"]
    print(f"pool depth: {pool_sizes.min()} / {pool_sizes.median():.0f} / {pool_sizes.max()} "
          f"(min / median / max)")
    return CA, RK, fetch_file, held, np, pl


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Scoring the pool

    The unit a cassette carries is the **long peptide around the mutation** — a 25-mer here — not a
    minimal epitope: a 9-mer loads onto any cell without costimulation and is the tolerising
    configuration. Each channel is maximised over the unit's 8–11mer windows, which is what a
    set-level reading of a unit asks: does it *contain* a window far from self, a window that is
    viral-like.

    `pres` and `occupancy` are supplied as `nan`. That is not a workaround — `aggregate_score`
    documents that a non-finite value takes the training mean, which is what the fit itself did for a
    candidate missing a column, and this table carries no HLA typing. The cost is measured: AUROC
    0.7246 CD4⁺ / 0.7216 CD8⁺ on this half **with no genotype**, against 0.7379 / 0.6907 for
    published baselines that had one.
    """)
    return


@app.cell
def _(CA, RK, held, np, pl):
    from mhcmatch.cli import _aggregate_channels
    from mhcmatch import complement as CM

    def _windows(wt, mt):
        """Every 8-11mer of `mt` spanning the first residue where it differs from `wt`."""
        if not mt:
            return []
        pos = next((i for i in range(min(len(wt or ""), len(mt))) if wt[i] != mt[i]), 0)
        return sorted({mt[s:s + L] for L in (8, 9, 10, 11)
                       for s in range(max(0, pos - L + 1), min(pos, len(mt) - L) + 1)
                       if len(mt[s:s + L]) == L})

    _wins = pl.DataFrame([{"mutation_id": m, "w": w}
                          for m, wt, mt in held.select("mutation_id", "wt_seq", "mt_seq").iter_rows()
                          for w in _windows(wt or "", mt or "")],
                         schema={"mutation_id": pl.Utf8, "w": pl.Utf8}).unique(maintain_order=True)
    _peps = sorted(set(_wins["w"].to_list()))
    print(f"{len(_peps):,} distinct windows over {_wins['mutation_id'].n_unique():,} mutations")

    _ch = pl.DataFrame({"w": _peps})
    for _col, _scale in RK.PHYS_COLUMNS.items():
        _ch = _ch.with_columns(pl.Series(_col, CM.burial(_peps, scale=_scale)))
    _got = _aggregate_channels("mhc1", no_self=False, species="human")(_peps)
    for _col in RK.CHANNEL_COLUMNS:
        _ch = _ch.with_columns(pl.Series(_col, _got[_col]))

    _cols = [c for c in _ch.columns if c != "w"]
    _agg = (_wins.join(_ch, on="w", how="left").group_by("mutation_id")
                 .agg(*[pl.col(c).max().alias(c) for c in _cols]))
    # `expr_pct` is a WITHIN-DONOR percentile -- the definition it has in the fit and in
    # `bench/cassette/select.py`, so the numbers below line up with `cassette_select.md`. The source
    # column is an RNA quartile 1-4; ranking it inside the donor spreads the ties onto (0, 1].
    scored = held.join(_agg, on="mutation_id", how="left").with_columns(
        pl.col("rna_exp_qrt").cast(pl.Float64, strict=False).alias("_ab")).with_columns(
        pl.when(pl.col("_ab").is_not_null())
          .then((pl.col("_ab").rank("average").over("patient") - 0.5)
                / pl.col("_ab").is_not_null().sum().over("patient"))
          .otherwise(0.5).alias("expr_pct"))
    feats = {f: (scored[f].to_numpy().astype(float) if f in scored.columns
                 else np.full(scored.height, np.nan)) for f in RK.AGGREGATE_FEATURES}
    scored = scored.with_columns(pl.Series("epic", RK.aggregate_score(feats))) \
                   .filter(pl.col("epic").is_finite()).sort(["patient", "mutation_id"])
    print(f"EPIC over {scored.height:,} units: "
          f"{scored['epic'].mean():+.3f} +/- {scored['epic'].std():.3f}")
    return (scored,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. One offset over the corpus, held

    `rank.probability` anchors the mean of **the batch it is handed**. Fitted once here, over every
    donor at once, so two donors' probabilities are on one axis. The alternative — one call per donor
    — is section 5.
    """)
    return


@app.cell
def _(CA, RK, scored):
    b = CA.prob_offset(scored["epic"].to_numpy().astype(float), RK.POOL_PREVALENCE)
    print(f"offset b = {b:+.4f} at prevalence {RK.POOL_PREVALENCE:.4f} "
          f"({RK.POOL_PREVALENCE * 615:.0f} of 615, TESLA)")
    return (b,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. `select` against a sort, on the objective each maximises

    Three arms on identical inputs, so a difference is the rule and nothing else. **`H` is built once
    over each donor's pool and every arm evaluated against it** — `goal_energy` renormalises the
    overlap to the set it is handed, so a per-cassette `H` would put each arm on its own basis and
    cancel exactly what a selection rule buys. `cassette.score` declines to report `H` for that
    reason; this is the five lines its docstring points at.
    """)
    return


@app.cell
def _(CA, b, np, pl, scored):
    K = 10
    _rng = np.random.default_rng(0)
    _rows = []
    for (_pat,), _g in scored.group_by("patient", maintain_order=True):
        if _g.height < K + 2:
            continue
        _s = _g["epic"].to_numpy().astype(float)
        _peps = _g["mt_seq"].to_list()
        _lab = ((_g["cd8"] + _g["cd4"]) > 0).to_numpy().astype(int)

        _c = CA.select(_s, _peps, None, k=K)
        _p = 1.0 / (1.0 + np.exp(-np.clip(_s + _c.offset, -60, 60)))
        _h, _J = CA.goal_energy(_p, CA.overlap(_peps, strength=_s), rho=_c.rho)
        _arms = {"select": _c.index,
                 "sort": list(np.argsort(-_s, kind="stable")[:K]),
                 "random": list(_rng.choice(_g.height, size=K, replace=False))}
        for _arm, _idx in _arms.items():
            _sc = CA.score(_s, _peps, None, chosen=list(_idx),
                           pool_scores=_s, pool_peptides=_peps, offset=b)
            _rows.append({"patient": _pat, "arm": _arm, "H": CA.energy(_h, _J, list(_idx)),
                          "yield": _sc["yield"], "lam": _sc["lam"],
                          "rho_dom": _sc["rho_dom"], "captured": int(_lab[list(_idx)].sum()),
                          "pool_pos": int(_lab.sum()), "pool_n": _g.height})
    arms_df = pl.DataFrame(_rows)
    summary = (arms_df.group_by("arm")
                      .agg(pl.len().alias("donors"),
                           pl.col("H").median().alias("H_median"),
                           pl.col("yield").median().alias("yield_median"),
                           pl.col("lam").median().alias("lam_median"),
                           pl.col("rho_dom").median().alias("rho_dom_median"),
                           pl.col("captured").sum().alias("captured"))
                      .sort("arm"))
    n_don = arms_df["patient"].n_unique()
    pos = int(arms_df.filter(pl.col("arm") == "select")["pool_pos"].sum())
    tot = int(arms_df.filter(pl.col("arm") == "select")["pool_n"].sum())
    print(summary)
    print(f"\n{n_don} donors, {tot:,} candidates, {pos} responding; base rate {pos / tot:.4f}")
    print(f"a rule that knew nothing would take {pos / tot * K * n_don:.1f} of them")
    return K, arms_df, n_don, pos, summary, tot


@app.cell
def _(mo):
    mo.md(r"""
    Read the two medians as the two different objectives they are. **The sort wins on `yield` by
    construction**, because that is what a sort maximises; **`select` wins on `H`**, which is what it
    maximises, and it pays for that in expected count. The `rho_dom` column is where the difference is
    spent: a larger mean strength gap between two units of the cassette, which is the immunodominance
    axis written as a pairwise quantity.

    Neither column is *the* answer. The corpus where the joint question can actually be settled is
    IVAC MUTANOME, where every one of 125 manufactured units was assayed — and there the two rules
    take 36 and 35 responding units of 60, agreeing on 45 of 60 slots.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 5. The calibration offset decides what is being reported

    One offset over the corpus, against one per donor, on the same cassettes.
    """)
    return


@app.cell
def _(CA, RK, K, np, pl, scored):
    _per = []
    for (_pat,), _g in scored.group_by("patient", maintain_order=True):
        if _g.height < K + 2:
            continue
        _s = _g["epic"].to_numpy().astype(float)
        _peps = _g["mt_seq"].to_list()
        _c = CA.select(_s, _peps, None, k=K)
        _shared = CA.score(_s, _peps, None, chosen=_c.index, offset=_c.offset)
        _own = CA.score(_s, _peps, None, chosen=_c.index, prevalence=RK.POOL_PREVALENCE)
        _mean = float(np.mean(1 / (1 + np.exp(-(_s + CA.prob_offset(_s, RK.POOL_PREVALENCE))))))
        _per.append({"patient": _pat, "pool_n": _g.height, "pool_mean_p_own": _mean,
                     "yield_batch": _shared["yield"], "yield_own": _own["yield"]})
    cal = pl.DataFrame(_per)
    print(f"pools span {cal['pool_n'].min()} to {cal['pool_n'].max()} candidates\n")
    print("one offset per donor -- every donor's POOL mean:")
    print(f"  min {cal['pool_mean_p_own'].min():.9f}   max {cal['pool_mean_p_own'].max():.9f}   "
          f"sd {cal['pool_mean_p_own'].std():.3e}")
    print(f"\ncassette `yield`, batch offset: {cal['yield_batch'].min():.3f} "
          f"to {cal['yield_batch'].max():.3f}  (spread {cal['yield_batch'].std():.4f})")
    print(f"cassette `yield`, own offset  : {cal['yield_own'].min():.3f} "
          f"to {cal['yield_own'].max():.3f}  (spread {cal['yield_own'].std():.4f})")
    return (cal,)


@app.cell
def _(mo):
    mo.md(r"""
    Every donor's pool mean lands on the declared prevalence to machine precision, whatever their
    pool holds. Read as a probability that number is not one, and two donors' numbers are the same
    number — which is why `MHCMATCH_CASSETTE_SCORE` in the Nextflow module collects every sample
    before scoring, and is the only process in that subworkflow that is deliberately not per sample.

    It is a real quantity, though, and the stronger of the two against immune infiltrate: on 4,073
    TCGA donors across 30 tumour types the per-donor-anchored sum reaches ρ = **+0.1298** against
    **+0.1115** for the batch-anchored one. It is an **enrichment** — how far a donor's chosen units
    sit above their own background — rather than a level. Ask for it deliberately:
    `cassette.group_offsets`, or `mhcmatch cassette score --per-donor-offset`.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6. What a trial actually manufactured

    `vaccines/cassettes.parquet` is the composition of every published NEO-PV-01 cassette — 1,506
    synthetic long peptides over 81 patients of two trials — plus IVAC MUTANOME's 125 RNA-pentatope
    units. No candidate pool was published for any of them, so `lam` is not computable there and
    `score` returns `None` for it rather than inventing a comparison.
    """)
    return


@app.cell
def _(CA, fetch_file, np, pl):
    _cass = pl.read_parquet(fetch_file("vaccines/cassettes.parquet"))
    _iv = pl.read_parquet(fetch_file("vaccines/ivac_mutanome_units.parquet"))
    made = pl.concat([_cass.select("trial", "patient", "peptide"),
                      _iv.select("trial", "patient", "peptide")], how="vertical")
    sz = made.group_by(["trial", "patient"]).len()
    print(made.group_by("trial").agg(pl.col("patient").n_unique().alias("cassettes"),
                                     pl.len().alias("units")).sort("trial"))
    print(f"\ncassette size across the pooled corpus: {sz['len'].min()} to {sz['len'].max()} units "
          f"-- a range no single trial contains")
    print(f"peptide length: {made['peptide'].str.len_chars().min()} to "
          f"{made['peptide'].str.len_chars().max()} residues")
    return made, sz


@app.cell
def _(mo):
    mo.md(r"""
    ## What to take away

    - **Give `select` the whole pool.** `expr_pct` (+0.3007) and `pres` (+0.2200) carry the two
      largest coefficients in the shipped model, so a shortlist already cut on binding and expression
      has no range left along them. Measured: on this exhaustive screen, responding at 0.0144,
      selecting five units captures 3.92× the base rate. On TESLA's *nominated* list, responding at
      0.0612 — 4.25× the same rate — every rule sits at the base rate, because the selection had
      already been done to it.
    - **Compare rules on `H`, not on `yield`.** A sort maximises `yield` and will always win on it.
    - **Compare cassettes on `lam`.** It divides by the donor's own pool, so it crosses donors and
      sizes without any shared calibration.
    - **Decide which offset you want.** A batch offset gives a level; a per-donor offset gives an
      enrichment. Both are useful and they are not the same number.
    """)
    return


if __name__ == "__main__":
    app.run()
