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
    # 12 — Allotype coverage, HLA loss, and tumour selectivity (`mhcmatch.cassette`)

    **What this demonstrates.** Three things a designer asks that a per-unit score cannot answer,
    on the six TESLA donors — 605 nominated candidates with 37 validated-immunogenic among them,
    the one public corpus that publishes a real per-candidate restriction *and* a measured label.

    1. **What allotype coverage is, and why the denominator is the donor's genotype.** An allotype
       holding **zero** units is the inequality the index exists to report, and it is exactly the
       one a coverage taken over the cassette's own labels cannot see.
    2. **What pricing HLA loss buys.** Every unit on one class-I molecule is lost together if that
       allele is. `block_live` puts that in the objective as the *exact* covariance it implies, and
       the readout is `captured_loh`: validated units still in the cassette after the worst single
       allotype is lost.
    3. **Tumour selectivity is a stated preference, not a fitted one.** The shipped model fits
       `expr_norm` — the source gene's level in *healthy* tissue — at **+0.4950**, its largest
       positive coefficient, because it answers *will this respond*. "High in tumour, low in normal"
       is a different, safety question, and it enters as a declared exchange rate.

    Everything is fetched from the public `isalgo/pmhc_data` deposit, so this runs on a bare
    `pip install mhcmatch`. Roughly 60 s, most of it the corpus channels.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. Six donors, their allotypes, and their answer key

    `neoantigens/neoag_tested.tsv.gz` carries TESLA as 605 rows with a per-candidate class-I
    restriction and `immunogenicity` as the measured label. That restriction is what makes this
    corpus able to answer an allotype question at all — the NCI screen, which has four times the
    donors, publishes no genotype.
    """)
    return


@app.cell
def _():
    import numpy as np
    import polars as pl

    from mhcmatch import cassette as CA
    from mhcmatch import rank as RK
    from mhcmatch.store import fetch_file

    _all = pl.read_csv(fetch_file("neoantigens/neoag_tested.tsv.gz"), separator="\t",
                       infer_schema_length=5000)
    tesla = (_all.filter(pl.col("dataset_origin") == "TESLA")
                 .select(donor=pl.col("patient_id").cast(pl.Utf8),
                         peptide=pl.col("peptide").str.to_uppercase(),
                         allele=pl.col("mhc_a").cast(pl.Utf8),
                         abundance=pl.col("expression").cast(pl.Float64, strict=False),
                         label=pl.col("immunogenicity").cast(pl.Int64))
                 .unique(subset=["donor", "peptide", "allele"], maintain_order=True)
                 .sort(["donor", "peptide"]))
    tesla.group_by("donor").agg(
        pl.len().alias("candidates"), pl.col("label").sum().alias("validated"),
        pl.col("allele").n_unique().alias("allotypes")).sort("donor")
    return CA, RK, fetch_file, np, pl, tesla


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Coverage, and the denominator that decides what it means

    `portfolio.coverage` returns per-allotype counts, a Gini index, the share of maximum entropy,
    and `n_covered` of `n_allotypes`. **`universe` is the donor's *distinct* allotypes**, and
    passing it is the whole point in both directions:

    * a patient homozygous at *B* has five distinct class-I allotypes, not six, so an even cassette
      over five is perfectly even — scoring it against six would report a genotype as a design flaw;
    * an allotype in `universe` with no unit is counted as a **zero**, which is the inequality the
      index is for.
    """)
    return


@app.cell
def _(pl):
    from mhcmatch import portfolio as PF

    _cass = ["A*02:01", "A*02:01", "A*02:01", "B*44:02"]
    _rows = [
        {"denominator": "the cassette's own labels", **{
            k: v for k, v in PF.coverage(_cass).items() if k != "counts"}},
        {"denominator": "the donor's genotype (6 allotypes)", **{
            k: v for k, v in PF.coverage(
                _cass, universe=["A*02:01", "A*03:01", "B*44:02", "B*07:02",
                                 "C*05:01", "C*07:01"]).items() if k != "counts"}},
        {"denominator": "a donor homozygous at B (5 allotypes)", **{
            k: v for k, v in PF.coverage(
                _cass, universe=["A*02:01", "A*03:01", "B*44:02",
                                 "C*05:01", "C*07:01"]).items() if k != "counts"}},
    ]
    pl.DataFrame(_rows)
    return (PF,)


@app.cell
def _(mo):
    mo.md(r"""
    The same four units read as **perfectly even** against their own labels (`entropy_ratio` is a
    share of the maximum over the allotypes present), and as **2 of 6 covered** against the
    genotype. Only the second is a statement about the cassette a patient receives.

    ## 3. Score the pool

    Genotype-free, the convention `bench/results/cassette_select.md` uses: `binder` and `log10a`
    are supplied as `nan`, so the shipped `aggregate_score` takes its documented training-mean path
    and the only difference between the arms below is the **selection rule**.
    """)
    return


@app.cell
def _(CA, RK, np, pl, tesla):
    from mhcmatch.cli import _aggregate_channels
    from mhcmatch.complement import burial

    # The corpus channels come from the builder `mhcmatch rank` itself scores with, never from a
    # local contraction: a channel computed against a different geometry is not a smaller effect,
    # it is a different feature, and `aggregate_score` says so in its own docstring.
    _peps = sorted(set(tesla["peptide"].to_list()))
    _got = _aggregate_channels("mhc1", no_self=False, species="human")(_peps)
    _feat = {c: dict(zip(_peps, burial(_peps, scale=sc))) for c, sc in RK.PHYS_COLUMNS.items()}
    _feat.update({c: dict(zip(_peps, _got[c])) for c in RK.CHANNEL_COLUMNS})
    scored = tesla.with_columns(
        [pl.col("peptide").replace_strict(m, default=None).alias(f) for f, m in _feat.items()])

    # `expr_lvl = log2(1 + TPM/c)`. TESLA publishes a real tumour TPM, so this term is defined here;
    # `c` is a property of a transcriptome and comes from the reference, never from the batch. The
    # pooled floor, because this table's `cancer_type` names the assayed compartment (PBMC / TIL)
    # rather than a tumour type. `expr_norm` needs a gene symbol, which TESLA does not publish, so
    # it takes `aggregate_score`'s documented training-mean path -- stated, not imputed.
    from mhcmatch.expression import context_floor

    _c = float(context_floor())
    scored = scored.with_columns(
        pl.when(pl.col("abundance").is_not_null())
          .then((1.0 + pl.col("abundance").clip(0.0) / _c).log(2))
          .alias("expr_lvl"))
    print(f"expr_lvl on {scored['expr_lvl'].is_not_null().sum()}/{scored.height} rows "
          f"at floor c = {_c:.4f} TPM; expr_norm on 0 (no gene symbol in this deposit)")
    _cols = {f: (scored[f].to_numpy().astype(float) if f in scored.columns
                 else np.full(scored.height, np.nan)) for f in RK.AGGREGATE_FEATURES}
    scored = scored.with_columns(pl.Series("epic", RK.aggregate_score(_cols)))
    _b = CA.prob_offset(scored["epic"].to_numpy().astype(float), RK.POOL_PREVALENCE)
    print(f"offset b = {_b:+.4f} at prevalence {RK.POOL_PREVALENCE:.4f} "
          f"over {scored.height} candidates")
    return (scored,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. Three rules, and what survives losing an allotype

    `sort` is the top *k* by score — the cassette a ranked list gives you. `select` is the shipped
    mean–variance objective. `select+loh` is the same objective with `block_live = 0.8`, which adds
    to `J_ij` exactly

    $$\mathrm{Cov}(R_i, R_j) = (1 - q)\, p_i p_j / q$$

    for two units on one allotype and nothing across allotypes — the covariance a lost allele
    implies, derived rather than fitted. **`q = 1.0` is bit-identical to no parameter at all.**

    `captured` is validated units in the cassette. `captured_loh` is how many are left after the
    **worst single** allotype is lost: the bad draw, not an average over draws, because loss of
    heterozygosity takes a specific allele.
    """)
    return


@app.cell
def _(CA, RK, np, pl, scored):
    def arms(g, k):
        s = g["epic"].to_numpy().astype(float)
        peps, alle = g["peptide"].to_list(), g["allele"].to_list()
        return {"sort": list(np.argsort(-s, kind="stable")[:k]),
                "select": CA.select(s, peps, alle, k=k).index,
                "select+loh": CA.select(s, peps, alle, k=k, block_live=0.8).index}

    _off = CA.prob_offset(scored["epic"].to_numpy().astype(float), RK.POOL_PREVALENCE)
    _rows = []
    for _k in (10, 20):
        for (_donor,), _g in scored.group_by("donor", maintain_order=True):
            _g = _g.filter(pl.col("epic").is_finite()).sort("peptide")
            _lab, _al = _g["label"].to_numpy(), np.array(_g["allele"].to_list())
            for _arm, _idx in arms(_g, _k).items():
                _idx = list(_idx)
                _sc = CA.score(_g["epic"].to_numpy().astype(float), _g["peptide"].to_list(),
                               _g["allele"].to_list(), chosen=_idx, offset=_off,
                               universe=sorted(set(_al)))
                _a = _al[_idx]
                _rows.append({
                    "k": _k, "arm": _arm, "donor": _donor, "captured": int(_lab[_idx].sum()),
                    "captured_loh": min(int(_lab[_idx][_a != b].sum()) for b in set(_a)),
                    "yield": _sc["yield"], "yield_loh": _sc["yield_loh"],
                    "rho_hla": _sc["rho_hla"], "n_covered": _sc["coverage"]["n_covered"]})
    per_donor = pl.DataFrame(_rows)
    per_donor.group_by(["k", "arm"]).agg(
        pl.col("captured").sum(), pl.col("captured_loh").sum(),
        pl.col("yield").mean().round(4), pl.col("yield_loh").mean().round(4),
        pl.col("rho_hla").mean().round(4)).sort(["k", "arm"])
    return (per_donor,)


@app.cell
def _(mo):
    mo.md(r"""
    Read `yield_loh` and `rho_hla` first — they are levels over all six donors and every unit, so
    they are the stable columns. `rho_hla`, the share of pairs sharing an allotype, falls from about
    0.46 for the sort to about 0.28 once the loss is priced, and `yield_loh` — expected responding
    units surviving the worst draw — rises by roughly half. `captured` and `captured_loh` are
    **counts over 37 validated units across six donors**, so single-unit differences between arms
    here are within noise; the recorded per-donor table with all three sizes and the coverage floor
    is `bench/results/cassette_tesla_donors.md`, and it scores through the benchmark's own feature
    path rather than this notebook's abbreviated one, so its counts are not expected to match here
    unit for unit.

    Note that `select` already spread without being told to: the allotype channel was always one of
    the three in `overlap`. What `block_live` changes is that the spread stops being a side effect of
    an averaged heuristic and becomes a stated design parameter with a rate on it — and that the
    coupling is the covariance rather than one third of a mean.

    ## 5. Per donor, which is the unit a designer actually has
    """)
    return


@app.cell
def _(per_donor, pl):
    per_donor.filter(pl.col("k") == 20).sort(["donor", "arm"])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6. Tumour selectivity is stated, not fitted

    The shipped EPIC artifact fits both expression terms **positive**:

    | term | what it is | coefficient |
    |---|---|--:|
    | `expr_lvl` | this candidate's own source-gene abundance | +0.3704 |
    | `expr_norm` | the same gene's median in the tumour's **matched normal** tissue | **+0.4950** |

    That is not a defect. It was fitted on *will this respond*, and a gene transcribed everywhere
    responds more often — the largest coefficient in the model is telling you something true about
    immunogenicity. It is simply not the safety question, and the two must not be conflated:
    imposing the tumour/normal ratio on the fit (equal and opposite coefficients) would assert an
    answer the data rejects.

    So selectivity enters as a **declared exchange rate**, the way `gamma` does:

    $$h_i = p_i - \tfrac{\gamma}{2} s_i^2 + w\,(\mathrm{expr\_lvl}_i - \mathrm{expr\_norm}_i)$$

    with `w` in expected responding units per **log2-fold**. Two properties matter:

    * it is charged to the **objective**, never to `p` — `p` is a calibrated marginal that
      `portfolio.survival` reads literally, so discounting it would silently restate the response
      model as well as the preference;
    * `select` reports the trade, so `w` is auditable rather than a knob.
    """)
    return


@app.cell
def _(CA, np):
    _rng = np.random.default_rng(12)
    _n = 120
    _s = _rng.normal(-1.0, 1.5, _n)
    _peps = ["".join(_rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 9)) for _ in range(_n)]
    _al = list(_rng.choice(["A*02:01", "B*44:02", "C*05:01"], _n))
    _lvl, _nrm = _rng.uniform(0, 8, _n), _rng.uniform(0, 8, _n)
    _d = CA.selectivity_delta(_lvl, _nrm)

    print(f"{'w':>6}  {'yield':>7}  {'mean log2 tumour/normal':>24}  slots changed")
    _base = CA.select(_s, _peps, _al, k=15)
    for _w in (0.0, 0.02, 0.05, 0.10):
        _c = CA.select(_s, _peps, _al, k=15, selectivity=_w, expr_lvl=_lvl, expr_norm=_nrm)
        print(f"{_w:>6.2f}  {_c.yield_:>7.4f}  {_d[_c.index].mean():>24.3f}  "
              f"{len(set(_c.index) - set(_base.index)):>2} of 15")
    return


@app.cell
def _(mo):
    mo.md(r"""
    Every reported `p` is identical across those rows — the same `sigma(s + b)` map — because the
    weight never touched it. What moved is which units the objective chose, and the run says what
    that cost in expected responding units.

    ## What to take away

    1. **Coverage without a `universe` answers a different question.** An allotype holding zero
       units is invisible against the cassette's own labels, and a homozygous donor is scored as a
       design flaw against a denominator of six.
    2. **An allotype is a group of units that fail together**, so it belongs in the objective, and
       the block response model says exactly how: `(1 - q) p_i p_j / q` on same-allotype pairs and
       zero elsewhere. Nothing is fitted; `q` is stated.
    3. **`captured_loh` and `yield_loh` are what to design against** when HLA loss is the worry.
       Two cassettes with the same `captured` and different `captured_loh` are different objects,
       and only the second column can tell them apart.
    4. **A fitted coefficient answers the question it was fitted on.** `expr_norm` at +0.4950 is a
       fact about immunogenicity, not permission to ignore tumour selectivity — and not a reason to
       refit either. State the preference, charge it to the objective, and report the trade.

    The per-donor benchmark behind §4, with the coverage floor and all three sizes, is
    `bench/results/cassette_tesla_donors.md` in the benchmark repository.
    """)
    return


if __name__ == "__main__":
    app.run()
