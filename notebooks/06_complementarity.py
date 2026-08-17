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
    # 6 — Complementarity, on a whole corpus (`mhcmatch.complement`)

    **What this demonstrates.** The recognition axis — *of the peptides an allele does present, which
    ones can a T-cell repertoire see* — as one score with six feature blocks, run over an entire
    published corpus rather than a hand-picked handful.

    **What you should conclude.** Two things, and the second is the surprising one:

    1. The blocks are not interchangeable descriptors. Where a residue sits (buried vs facing the
       receptor), what kind of statistic is taken (a property average vs a contiguous motif vs
       residue identity), and which contact potential applies to which side are all different
       questions, and the answers disagree.
    2. **Residue identity carries almost all of it.** The physicochemical blocks — the ones
       `ipred` is built from — are real but small next to a per-role log-odds over the 20 amino
       acids. The notebook shows this rather than asserting it.

    Everything here is vendored except the corpus in section 6.4.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""## 6.1 — One peptide, and its parts""")
    return


@app.cell
def _():
    from mhcmatch import complement

    peps = ["GILGFVFTL", "NLVPMVATV", "SIINFEKL", "KRWIILGLNK", "AAAAAAAAA"]
    scores = complement.score(peps)
    dict(zip(peps, [round(float(s), 4) for s in scores]))
    return complement, peps, scores


@app.cell
def _(complement):
    # The design matrix, one row per peptide, in the order the coefficients expect.
    complement.features("GILGFVFTL")
    return


@app.cell
def _(mo):
    mo.md(r"""
    The score is a **log-odds and carries no prior** — the same contract as
    `mhcmatch.posbayes.llr`. The training corpus runs at ~3.2% positives, a viral proteome scan
    nearer 3.0e-3, the NCI neoantigen screen 4.2e-4. Reading a corpus-prevalence probability as an
    operational one overstates it by up to 75x, so the base rate is yours to supply.
    """)
    return


@app.cell
def _(complement, peps):
    corpus_rate = complement.posterior(peps, prior=complement.PARAMS["prevalence"])
    screen_rate = complement.posterior(peps, prior=4.2e-4)
    {p: (round(float(a), 4), float(f"{b:.3g}")) for p, a, b in zip(peps, corpus_rate, screen_rate)}
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6.2 — Arrangement is not composition

    A run of 3–4 hydrophobic residues facing the TCR is a different object from the same residues
    scattered along the peptide, and **no sum can express the difference**: a permutation has
    identical sums by construction. The `motif` block exists for exactly this.

    TCR-facing positions of a 9-mer are 3..6 under the class-I `pockets` split, so `IIDD` is one run
    of 2 and `IDID` is two runs of 1 — same composition, same `run_frac`, different arrangement.
    """)
    return


@app.cell
def _(complement):
    f_arr, _ = complement.encode(["AAAIIDDAA", "AAAIDIDAA"])
    {k: f_arr[k].tolist() for k in ("kd_run_max", "kd_run_n", "kd_run_frac", "pc1_tcr")}
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6.3 — The `aa` block *is* `posbayes`

    The per-role residue log-odds columns sum to exactly `mhcmatch.posbayes.llr`. That is not a
    coincidence to be checked at review time — it means the shipped position-role model is this
    feature set's `aa` block on its own, so the block ablation measures what everything else adds
    to a model that already exists.
    """)
    return


@app.cell
def _(complement, peps):
    import numpy as np

    from mhcmatch import posbayes

    _, ct = complement.encode(peps)
    t = posbayes.table("human")
    mine = ct["anchor"] @ np.array(t["anchor"]) + ct["tcr"] @ np.array(t["tcrface"])
    {p: (round(float(m), 6), round(posbayes.llr(p), 6)) for p, m in zip(peps, mine)}
    return (np,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 6.4 — A whole corpus, in one call

    `score` is vectorised: the whole feature set is two `(n, 20)` count matrices times a handful of
    property vectors, plus a sparse pair list. Hand it everything at once — looping peptide by
    peptide is the slow path and buys nothing.

    `store.fetch_file` pulls any file of the public `isalgo/pmhc_data` dataset and caches it. Set
    `MHCMATCH_PMHC_DIR` to a local mirror to skip the download.
    """)
    return


@app.cell
def _(complement, np):
    import csv
    import gzip
    import time

    from mhcmatch import store

    path = store.fetch_file("immunogenicity/chowell_rebuilt.tsv.gz")
    with gzip.open(path, "rt") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    corpus_peps = [r["peptide"] for r in rows]
    y = np.array([int(r["label"]) for r in rows])
    t0 = time.time()
    s = complement.score(corpus_peps)
    elapsed = time.time() - t0
    (len(rows), int(y.sum()), round(elapsed, 2))
    return corpus_peps, s, y


@app.cell
def _(mo):
    mo.md(r"""
    Immunogenic peptides score higher, and the separation is worth quantifying rather than eyeballing.
    This is an **in-sample** figure — this corpus is what the shipped parameters were fitted on — so
    it is a sanity check, not a performance claim. The out-of-sample numbers (peptide-grouped
    cross-validation, corpus transfer, cross-species transfer) live in the benchmark repository at
    `bench/results/complementarity.md`.
    """)
    return


@app.cell
def _(np, s, y):
    def auroc(pos, neg):
        r = np.argsort(np.argsort(np.concatenate([pos, neg]))) + 1
        return (r[: len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))

    {
        "mean score, immunogenic": round(float(s[y == 1].mean()), 4),
        "mean score, presented-not-immunogenic": round(float(s[y == 0].mean()), 4),
        "AUROC (IN-SAMPLE)": round(float(auroc(s[y == 1], s[y == 0])), 4),
    }
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 6.5 — From the command line

    ```bash
    mhcmatch complement GILGFVFTL SIINFEKL NLVPMVATV

    # the whole deposit; --peptides takes one-per-line or a TSV with a `peptide` column
    mhcmatch complement --peptides chowell_rebuilt.tsv.gz --prior 3.2e-2 --out scored.tsv

    # every feature, so a score can be taken apart
    mhcmatch complement GILGFVFTL --features
    ```

    ## Caveat

    **Class I only.** The role split is the class-I one (P1–P3, PΩ-1, PΩ). A class-II ligand is
    anchored by the P1/P4/P6/P9 core of a 9-mer register floating inside a longer peptide, so
    applying this scheme to it labels the wrong residues as anchors and returns a confident, wrong
    number. `mhcmatch.rank` returns `NaN` for class II rather than guessing.
    """)
    return


if __name__ == "__main__":
    app.run()
