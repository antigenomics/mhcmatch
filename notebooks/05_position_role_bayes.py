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
    # 5 — Position-role immunogenicity (`mhcmatch.posbayes`)

    **What this demonstrates.** `ipred` scores a peptide from pooled physicochemical descriptors and
    does not distinguish *where* a residue sits. But an anchor residue is buried in the MHC groove and
    a TCR-facing residue is contacted by the receptor — different channels, and for several amino
    acids they carry **opposite signs**. Pooling averages that away.

    `posbayes` keeps them apart in the simplest form that can: for each role separately, the
    conditional amino-acid distribution given the class, scored as a summed log-likelihood ratio.

    Three things in this notebook are worth more than the score itself:

    1. the model emits a **log-likelihood ratio**, so you supply the prior — and the right prior is
       nowhere near the training corpus's;
    2. **cysteine is masked**, because the negatives are mass-spec eluted and Cys is under-detected by
       detection chemistry rather than biology;
    3. the two role tables genuinely disagree, which is the model's entire content.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## The measured performance

    Peptide-grouped 5-fold cross-validation — no peptide in both train and test — on the IEDB
    positive-T-cell-assay vs self-eluted-ligand corpus.

    | | human | mouse |
    |---|--:|--:|
    | rows | 464,310 | 47,203 |
    | immunogenic | 14,712 | 5,154 |
    | **AUROC (`posbayes`)** | **0.712** | **0.758** |
    | AUROC (`ipred`, *in-sample*) | 0.607 | 0.668 |

    Size-matched cross-species transfer, mean over 10 matched subsamples:

    | direction | AUROC | sd |
    |---|--:|--:|
    | human → mouse | **0.731** | 0.003 |
    | mouse → human | **0.692** | 0.000 |

    > `ipred`'s figures are **in-sample** — that corpus is its training set — so this is not a
    > like-for-like contest. It is quoted because an in-sample baseline that still loses is the
    > conservative direction, not because it is fair.
    """)
    return


@app.cell
def _():
    from mhcmatch import posbayes

    peptides = ["GILGFVFTL", "NLVPMVATV", "KLGGALQAK", "SLYNTVATL", "AAAAAAAAA", "KKKKKKKKK"]
    scores = {p: posbayes.llr(p) for p in peptides}
    scores
    return peptides, posbayes, scores


@app.cell
def _(mo):
    mo.md(r"""
    ## The prior is yours to supply, and it matters by 1–2 orders of magnitude

    `llr()` carries **no** prior. `posterior()` takes one and has no default, deliberately: the
    training corpus runs at ~3.2% positives, a viral proteome scan at ~3.0e-3 (counted from distinct
    9-mers against known epitopes), and the NCI screen at 4.8e-4. A default would silently impose one
    setting's base rate on every caller.

    Recalibration is exact and additive — `logit P = llr + log(prior/(1-prior))` — so it never
    changes the ranking, only what the number means.
    """)
    return


@app.cell
def _(peptides, posbayes):
    PRIORS = {"training corpus": 0.0317, "viral proteome scan": 3.0e-3, "NCI screen": 4.81e-4}
    table = {
        p: {name: round(posbayes.posterior(p, pr), 5) for name, pr in PRIORS.items()}
        for p in peptides
    }
    table
    return PRIORS, table


@app.cell
def _(mo):
    mo.md(r"""
    ## Why cysteine is masked

    The negatives in the training corpus are **mass-spectrometry-eluted** ligands, and cysteine is
    systematically under-detected in immunopeptidomics unless alkylated. Measured on that corpus,
    Cys-containing peptides are:

    | evidence source | Cys-containing |
    |---|--:|
    | `iedb_tcell` (assayed positives) | **11.59%** |
    | `iedb_mhc_ligand` (MS-eluted negatives) | 1.73% |
    | `thymus_ms` (MS-eluted negatives) | **0.17%** |

    A 6.5× depletion driven by platform, not biology. Fitted freely, Cys took the **single largest
    coefficient** in the model (+1.84 anchor / +2.05 TCR-facing) — the model was reading assay
    platform. Zeroing it costs grouped-CV AUROC 0.712 → 0.690 and *improves* transfer to neoantigen
    screens (TESLA 0.393 → 0.435, NCI 0.489 → 0.521).

    Any model trained on MS-eluted negatives against assayed positives inherits this, `ipred`
    included.
    """)
    return


@app.cell
def _(posbayes):
    c = posbayes.AA.index("C")
    cys_is_zero = {
        "human anchor": posbayes.HUMAN["anchor"][c],
        "human tcrface": posbayes.HUMAN["tcrface"][c],
        "mouse anchor": posbayes.MOUSE["anchor"][c],
        "mouse tcrface": posbayes.MOUSE["tcrface"][c],
    }
    cys_is_zero
    return c, cys_is_zero


@app.cell
def _(mo):
    mo.md(r"""
    ## The two role tables disagree — that is the whole point

    If anchor and TCR-facing agreed everywhere, splitting them would buy nothing. Residues where the
    sign flips between roles are the ones a pooled model must average into noise.
    """)
    return


@app.cell
def _(posbayes):
    flips = {
        aa: (round(posbayes.HUMAN["anchor"][i], 3), round(posbayes.HUMAN["tcrface"][i], 3))
        for i, aa in enumerate(posbayes.AA)
        if posbayes.HUMAN["anchor"][i] * posbayes.HUMAN["tcrface"][i] < 0
    }
    flips
    return (flips,)


@app.cell
def _(mo):
    mo.md(r"""
    ## Where this does *not* work, stated plainly

    On the neoantigen screens the model is weak to inverted — TESLA 0.435, NCI 0.521, Gfeller 0.656,
    Gfeller-GBM 0.550 — against presentation (`binder`) at 0.785 / 0.967 / 0.794 / 0.647.

    Two measured facts constrain the explanation, and neither has been resolved:

    * the inversion sits in the **anchor** half, and is not uniform (anchor-only is 0.403 on TESLA but
      0.620 on Gfeller — the same table pointing opposite ways on two cohorts);
    * a single amino-acid substitution moves the score by a median of 0.283 against a between-peptide
      spread of sd 1.096, so the model is largely blind to the mutation that *makes* a neoantigen.

    Use it as a recognition term on presented peptides, not as a standalone neoantigen ranker.
    """)
    return


if __name__ == "__main__":
    app.run()
