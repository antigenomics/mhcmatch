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
    # 9 — Composition is not ranking (`mhcmatch.portfolio`)

    **What this demonstrates.** Why a cassette is a *set* problem rather than a ranking problem, on
    synthetic data small enough that every number is checkable by hand. Notebook 8 builds a cassette;
    this one asks what the cassette is *worth*, and shows two things a scalar score cannot do.

    **What you should conclude.**

    1. Two cassettes with the **identical** expected number of responders can differ twofold in the
       probability that *at least one* works. The scalar objective cannot tell them apart.
    2. A weighted sum can only ever rank first what sits on the **upper convex hull** of the
       objective cloud. Some Pareto-efficient candidates are reachable by no weighting at all.
    3. Swapping in a Chebyshev aggregator fixes (2). Nothing fixes (1) short of changing the
       selection rule, because `P(>= k | S)` is not a sum over candidates.

    Nothing here is fitted. `mhcmatch.portfolio` takes the scores the rest of the library produces
    and reports what a proposed set of them is worth.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. The same expected yield, two different cassettes

    Twelve units, each with a marginal response probability of 0.15, so `sum(p) = 1.8` either way.
    The only difference is how they are spread over blocks — an allotype, a mechanism, anything that
    can fail as a unit for this donor. `q = 0.5` is the probability a block is live.
    """)
    return


@app.cell
def _():
    import numpy as np

    from mhcmatch import portfolio

    p = np.full(12, 0.15)
    one_block = np.zeros(12, dtype=int)          # all twelve on one allotype
    four_blocks = np.arange(12) % 4              # three each on four allotypes

    rows = []
    for name, blk in (("1 block", one_block), ("4 blocks", four_blocks)):
        n_b = int(blk.max()) + 1
        ge1 = portfolio.p_at_least(p, blk, [0.5] * n_b, k=1)
        ge2 = portfolio.p_at_least(p, blk, [0.5] * n_b, k=2)
        rows.append((name, float(p.sum()), ge1, ge2, portfolio.n_effective(p, ge1)))

    for name, sp, a, b, neff in rows:
        print(f"{name:>9}   sum(p)={sp:.2f}   P(>=1)={a:.3f}   P(>=2)={b:.3f}   n_eff={neff:.2f}")
    return four_blocks, np, one_block, p, portfolio


@app.cell
def _(mo):
    mo.md(r"""
    `sum(p)` is identical, `P(>=1)` is not, and `n_eff` says how many *independent* units each
    cassette is really worth. The one-block cassette is capped at `q = 0.5` however many units you
    add to it — that is the saturation proposition, and it is why adding a 31st strong unit on an
    allotype that already has thirty buys nothing.

    Check the cap directly:
    """)
    return


@app.cell
def _(np, portfolio):
    for m in (5, 20, 80):
        v = portfolio.p_at_least(np.full(m, 0.15), np.zeros(m, dtype=int), [0.5])
        print(f"m={m:>3}  all on one block  ->  P(>=1) = {v:.4f}   (cap q = 0.5)")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Selecting against blocks

    `vector.select` saturates a budget per block: it keeps adding on block `b` while
    `p_next > S_b / (n0 + n_b)`. Diversification falls out of that arithmetic rather than a quota.
    The default block is the allotype; `block=` generalises it. Here five candidates, two allotypes.
    """)
    return


@app.cell
def _():
    from mhcmatch.vector import Unit, select

    units = [
        Unit(peptide="AAAAAAAAA", mutation_index=4, gene="G0", allele="A*02:01", p=0.90, cls="mhc1"),
        Unit(peptide="CCCCCCCCC", mutation_index=4, gene="G1", allele="A*02:01", p=0.50, cls="mhc1"),
        Unit(peptide="DDDDDDDDD", mutation_index=4, gene="G2", allele="A*02:01", p=0.20, cls="mhc1"),
        Unit(peptide="EEEEEEEEE", mutation_index=4, gene="G3", allele="B*07:02", p=0.80, cls="mhc1"),
        Unit(peptide="FFFFFFFFF", mutation_index=4, gene="G4", allele="B*07:02", p=0.30, cls="mhc1"),
    ]
    corner_of = {"G0": "presentation", "G1": "recognition", "G2": "presentation",
                 "G3": "presentation", "G4": "recognition"}

    by_allele = select(units, n0=5.0)
    by_pair = select(units, n0=5.0, block=lambda u: (u.allele, corner_of[u.gene]))

    print("blocked on allotype      :", [u.gene for u in by_allele.units],
          f" expected_yield={by_allele.expected_yield:.3f}")
    print("blocked on allotype+corner:", [u.gene for u in by_pair.units],
          f" expected_yield={by_pair.expected_yield:.3f}")
    print("\nper block, allotype+corner:")
    for key, (n, s, y) in by_pair.per_block().items():
        print(f"  {str(key):<32} n={n}  sum(p)={s:.2f}  saturated={y:.3f}")
    return (select,)


@app.cell
def _(mo):
    mo.md(r"""
    `expected_yield` is computed against the partition the rule actually used — not
    unconditionally by allotype — or it would describe a cassette that was never built.

    `n0` is per-block capacity and has **no default**, deliberately: the dose-matched trial that
    would measure it does not exist in the public record. Sweeping it retrospectively on 178
    validated-immunogenic neoantigens puts the selection-layer optimum near 20.

    ## 3. What a weighted sum cannot select

    Five candidates, two objectives, higher is better on both. Candidate 4 is Pareto-efficient —
    nothing dominates it — yet it sits strictly inside the convex hull of candidates 0 and 2.
    """)
    return


@app.cell
def _(np, portfolio):
    Z = np.array([[3.0, 0.0],       # 0: extreme on objective 1
                  [0.0, 3.0],       # 1: extreme on objective 2
                  [1.6, 1.6],       # 2: balanced
                  [1.0, 1.0],       # 3: dominated by 2
                  [2.5, 0.5]])      # 4: efficient, but inside the hull

    front = portfolio.pareto_front(Z)
    print("Pareto-efficient   :", front)
    print("linearly supported :", [portfolio.linearly_supported(Z, i) for i in range(len(Z))])
    return (Z,)


@app.cell
def _(mo):
    mo.md(r"""
    Candidate 4 is non-dominated and yet supported by **no** `beta >= 0`. That is not a tuning
    problem: it has no supporting hyperplane, so no amount of weight search reaches it. On real
    candidate pools this is 45 of 161 Pareto-efficient validated neoantigens.

    A brute-force sweep confirms it, and shows the fix:
    """)
    return


@app.cell
def _(Z, np, portfolio):
    best_linear = np.full(len(Z), len(Z))
    for w0 in np.linspace(0, 1, 2001):
        scal = Z @ np.array([w0, 1 - w0])
        rank_w = 1 + (scal > scal[:, None]).sum(1)
        best_linear = np.minimum(best_linear, rank_w)
    print("best rank over 2001 weighted sums:", best_linear)

    witness = 4
    shortfall = (Z.max(0) + 1e-6) - Z[witness]
    lam = (1.0 / shortfall) / (1.0 / shortfall).sum()     # the closed-form Chebyshev witness
    cheb = portfolio.chebyshev_score(Z, lam)
    print(f"\nChebyshev with lambda ∝ 1/(z* - z_4) = {lam.round(3)}  ->  argmax = {cheb.argmax()}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    Rank 1 is unreachable for candidate 4 under every weighted sum, and immediate under Chebyshev,
    whose optimal weights are closed-form: `lambda_k ∝ 1 / (z*_k - z_k)`, which equalises the
    weighted shortfalls.

    **The part that does not get fixed.** A gradient-boosted score is not hull-limited either — an
    additive tree ensemble is dense in the continuous functions — and on real pools it is worth much
    more than any of this. But top-`m` by *any* pointwise score maximises `sum_i s_i`, a **modular**
    set function, while `P(>= k | S)` is submodular the moment two units share a block. That is a
    property of the selection *rule*, not the *scorer*, so no model capacity addresses it.

    ## 4. Measuring the dependence on your own readout

    Do not assume a correlation — measure it. One `(m, k)` per patient: units assayed, units
    positive. These are the counts from the adjuvant TNBC mRNA vaccine trial (Sahin et al.,
    *Nature* 2026;651:1088–1096).
    """)
    return


@app.cell
def _(portfolio):
    # the trial's own per-patient counts, unlabelled: 13 patients with an ex vivo bulk-PBMC
    # readout, most on the full 20-unit cassette. Two returned nothing; two returned eight.
    assayed = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 10, 5, 1]
    positive = [8, 8, 6, 3, 2, 2, 2, 2, 2, 1, 0, 5, 0]

    disp = portfolio.dispersion(assayed, positive)
    print(f"pooled per-unit rate   {disp['p_pooled']:.4f}")
    print(f"variance ratio         {disp['ratio']:.2f}x the independent-Bernoulli expectation")

    lrt = portfolio.betabinom_rho(assayed, positive)
    print(f"\nintra-patient rho      {lrt['rho']:.3f}")
    print(f"LRT vs binomial        D = {lrt['D']:.1f}, p = {lrt['p_value']:.2e}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    Keep the zero-response patients: they carry most of the information about dispersion, and a
    minimum-pool-size filter deletes precisely them. The p-value is conservative at these cohort
    sizes — simulation under the null puts the realised type-I error at 0.022 against a nominal 0.05.

    ### Where the numbers come from

    Every measured figure quoted here is generated in
    [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark) under
    `bench/vector/portfolio_*.py` and recorded in `bench/results/vector_portfolio_*.md`. The
    derivations are in the appendix chapter *Cassette composition as a portfolio*.
    """)
    return


if __name__ == "__main__":
    app.run()
