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
    # 7 — Mimicry as immune-response risk (`mhcmatch.mimicry`)

    **What this demonstrates.** Notebook 4 showed the raw scan. This one shows the *fitted* form:
    three references, each split into two channels, returning signed log-odds rather than counts.

    | component | reference | what a hit argues |
    |---|---|---|
    | `viral` | foreign presented ligandome | a pre-existing anti-pathogen repertoire may cross-react — **raises** expected immunogenicity |
    | `self` | the host proteome | tolerance, and simultaneously the **autoimmunity** flag |
    | `thymus` | thymic self-immunopeptidome | reactive precursors met it during negative selection |

    **What you should conclude.**

    1. A single whole-peptide distance is the wrong feature. Each component is split into an
       **anchor** channel and a **TCR-facing** channel, which partition the peptide.
    2. **Scores are log-odds. Probabilities need a named corpus.** The seven screens behind the
       shipped calibration run from 0.048 % to 46.8 % positive, so an unqualified
       `P(immunogenic)` is largely a statement about which intercept was used. Rank on log-odds.
    3. **The tested-neoantigen database is an annotation, never a fitted term** — §7.4. Every
       labelled screen we hold is inside it, so a coefficient on it would be memorisation.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 7.1 The shipped model, and what it says
    """)
    return


@app.cell
def _():
    from mhcmatch import mimicry

    p = mimicry.params("mhc1")
    print(f"model version {p['version']}   radius {p['radius']}")
    print(f"fitted on {p['fit']['n']:,} rows / {p['fit']['pos']:,} positives, "
          f"{len(p['fit']['screens'])} screens, screen indicators = {p['fit']['screen_indicators']}")
    print()
    for f, c, s in zip(p["features"], p["logistic"]["coef"], p["logistic"]["sd"]):
        print(f"  {f:<16}{c:+.4f}   z = {c / s:+6.2f}")
    return mimicry, p


@app.cell
def _(mo):
    mo.md(r"""
    Read the signs by **reference**: `viral` positive on both channels (priming), `self` negative on
    both (tolerance) — which is what the design predicted — and `thymus` positive on its anchor
    channel with the TCR channel unresolved.

    A different pattern (anchor positive, TCR-facing negative, across *every* reference) appears when
    mimicry is fitted **residual to a model that already contains a whole-peptide
    physicochemical term and a foreignness term**.
    That is a statement about what mimicry adds to those terms, not about mimicry alone. The two are
    easy to confuse and the module docstring keeps them apart deliberately.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 7.2 The two channels partition the peptide
    """)
    return


@app.cell
def _(mimicry):
    for L in (9, 10):
        m = mimicry.masks(L)
        print(f"L={L}  anchor {m['anchor']}   tcr {m['tcr']}")
        assert sorted(m["anchor"] + m["tcr"]) == list(range(L))
    print("\nno position is counted twice, and none is missed")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 7.3 Honest reporting: pooled vs within-screen

    The aggregate's own fit record carries both numbers, and the gap between them *is* the finding —
    pooling seven screens with different prevalences and candidate-generation processes manufactures
    AUROC.
    """)
    return


@app.cell
def _(p):
    print(f"AUROC pooled        {p['fit']['auroc_pooled']:.3f}")
    print(f"AUROC within screen {p['fit']['auroc_within_screen_median']:.3f}   <- report this one")
    print(f"\nblock LRT over screen indicators alone: chi2 = {p['fit']['lrt_chi2']:.2f} "
          f"on {p['fit']['lrt_df']} df, p = {p['fit']['lrt_p']:.3g}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 7.4 Prior evidence: the tested-neoantigen database

    `annotate` is a *lookup*, not a prediction, and it is kept out of the fitted aggregate on
    purpose. Use the fuzzy distance rather than exact matching: held out properly, matching at two
    substitutions recovers 0.08–0.34 of a fresh screen's positives where exact lookup recovers
    0.00–0.26.
    """)
    return


@app.cell
def _(mimicry):
    for r in mimicry.annotate(["KLVVVGACGV", "GILGFVFTL", "AAAWYLWEV"]):
        print(f"  {r['peptide']:<12} d={r['neoag_distance']}  known={r['known']!s:<6}"
              f"nearest={r['neoag_nearest'] or '-':<12} n_within={r['neoag_n_within']}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    Same thing from the shell, carrying every column of an existing candidate table through:

    ```
    mhcmatch neoag --peptides candidates.tsv --out candidates.annotated.tsv
    ```
    """)
    return


if __name__ == "__main__":
    app.run()
