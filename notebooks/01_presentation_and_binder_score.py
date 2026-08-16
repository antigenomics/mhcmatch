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
    # 1 — Presentation, restriction and the generalized binder score

    **What this demonstrates.** The core `mhcmatch` workflow end to end: build a `Store` from the
    public reference panel, split a peptide into its anchor / TCR-facing parts, rank the alleles that
    present it, and then score it with `binder_score` — the calibrated combination of the presentation
    head and the affinity head.

    **What you should conclude.** A `BinderScore` carries two different kinds of number and they answer
    different questions:

    * `binder_rank` is a **%rank** — where this (peptide, allele) pair sits against a random-peptide
      background for *that allele*. It is what you sort by.
    * `p_binder` is a **calibrated probability** on an absolute scale. It is what you threshold, hand to
      a downstream model, or compare across runs.

    The last two sections make the difference concrete: the *ordinal position* of a peptide in a
    candidate list moves as soon as the list changes, while its `p_binder` does not.

    Everything here bootstraps from the public HuggingFace dataset `isalgo/pmhc_data` — no local data
    files. First run downloads ~4 MB and caches it; later runs are instant.
    """)
    return


@app.cell
def _():
    import time

    import mhcmatch

    _t0 = time.time()
    # Auto-fetches pmhc/pmhc_shortlist.tsv.gz from the public HF dataset (cached by huggingface_hub).
    store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human")
    mhc1_alleles = store.alleles("mhc1")
    print(f"panel loaded in {time.time() - _t0:.1f} s")
    print(f"MHC-I human alleles in the shortlist tier: {len(mhc1_alleles)}")
    print(f"first five: {mhc1_alleles[:5]}")
    return mhc1_alleles, mhcmatch, store, time


@app.cell
def _(mo):
    mo.md(r"""
    ## 1.1 Anchor / TCR-facing decomposition

    `decompose` needs no reference data at all — it is pure layout. `X` marks the masked half, so the
    two strings read as "what the groove sees" and "what the TCR sees".
    """)
    return


@app.cell
def _(store):
    PEPTIDE = "NLVPMVATV"  # CMV pp65 495-503, the canonical HLA-A*02:01 epitope

    _d = store.decompose(PEPTIDE, cls="mhc1")
    print(f"peptide      {_d.peptide}")
    print(f"tcr_facing   {_d.tcr_facing}    (anchors masked -> the recognition read-out)")
    print(f"presentation {_d.presentation}    (TCR-facing masked -> the anchor read-out)")
    print(f"anchors      {_d.anchors}  (0-based)")
    return (PEPTIDE,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 1.2 Restriction — which allele presents this peptide?

    `restriction(calibrated=True)` ranks alleles on the **allele-specificity** axis (*which* allele,
    given that the peptide is presented at all) and attaches a per-allele `%rank`, a calibrated
    `p_present`, and a binding band.

    This is the first place a %rank and a probability appear side by side, and they already disagree in
    ordering — the %rank is computed against a random-peptide background, `p_present` against the
    allele's own known ligands.
    """)
    return


@app.cell
def _(PEPTIDE, mo, store, time):
    _t0 = time.time()
    restr = store.restriction(PEPTIDE, cls="mhc1", calibrated=True)
    print(f"restriction(calibrated=True) over the whole panel: {time.time() - _t0:.1f} s")


    def esc(s):
        """Escape the `*` in an HLA name so markdown does not read it as emphasis."""
        return s.replace("*", r"\*")


    _rows = "\n".join(
        f"| {esc(r.allele)} | {r.rank:.2f} | {r.p_present:.4f} | {r.band} | {r.n_votes} |"
        for r in restr[:6]
    )
    mo.md(
        f"""
    **Top presenting alleles for `{PEPTIDE}`**

    | allele | %rank | p_present | band | n_votes |
    |---|--:|--:|---|--:|
    {_rows}
    """
    )
    return (esc,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 1.3 `binder_score` — the generalized binder score

    `store.binder_score(peptide, alleles=..., cls=...)` returns a `BinderScore` per allele, sorted
    best-first. It runs two heads and fuses them:

    | field | head | meaning |
    |---|---|---|
    | `presentation_rank` | `AnchorModel` | %rank of the presentation log-odds vs a random-peptide background. Lower = stronger. |
    | `affinity_nm` | `PottsAffinity` | predicted IC50 in nM. Lower = stronger. |
    | `affinity_rank` | `PottsAffinity` | %rank of the affinity score for this allele. Lower = stronger. |
    | `binder_rank` | both | **calibrated combined %rank.** Fisher's statistic `-(ln p_pres + ln p_aff)` re-calibrated against the same random-peptide background, so the result is itself a true %rank. Lower = stronger. |
    | `band` | both | `strong` / `weak` / `non-binder`, banded on `binder_rank`. |
    | `p_binder` | both | **isotonic-calibrated `P(binder)`** over that same combined statistic, fit from the allele's known ligands against the background. Higher = stronger. |

    `p_binder` is populated because the internal calibrator is constructed with `positives=` (the
    allele's ligands from the panel), which is what enables the isotonic fit. Without positives the
    field would fall back to a rank-derived value.
    """)
    return


@app.cell
def _(PEPTIDE, esc, mo, store, time):
    SMALL_POOL = ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*01:01"]

    _t0 = time.time()
    small = store.binder_score(PEPTIDE, alleles=",".join(SMALL_POOL), cls="mhc1")
    print(f"binder_score over {len(SMALL_POOL)} alleles: {time.time() - _t0:.1f} s (cold calibrator)")

    _rows = "\n".join(
        f"| {i} | {esc(b.allele)} | {b.presentation_rank:.2f} | {b.affinity_nm:.1f} | "
        f"{b.affinity_rank:.2f} | {b.binder_rank:.2f} | {b.band} | {b.p_binder:.4f} |"
        for i, b in enumerate(small, 1)
    )
    mo.md(
        f"""
    **`{PEPTIDE}` against a 3-allele pool**

    | position | allele | presentation_rank | affinity_nm | affinity_rank | binder_rank | band | p_binder |
    |--:|---|--:|--:|--:|--:|---|--:|
    {_rows}
    """
    )
    return SMALL_POOL, small


@app.cell
def _(mo):
    mo.md(r"""
    ## 1.4 Widening the allele pool

    Now score the **same peptide** against a 12-allele pool that contains the three above. Two things
    to watch:

    * the per-allele `binder_rank` and `p_binder` of the shared alleles;
    * where `HLA-A*02:01` sits in the returned list.
    """)
    return


@app.cell
def _(PEPTIDE, SMALL_POOL, esc, mo, store, time):
    LARGE_POOL = SMALL_POOL + [
        "HLA-A*02:03",
        "HLA-A*02:06",
        "HLA-A*03:01",
        "HLA-A*11:01",
        "HLA-A*24:02",
        "HLA-B*08:01",
        "HLA-B*35:01",
        "HLA-B*44:02",
        "HLA-C*07:01",
    ]

    _t0 = time.time()
    large = store.binder_score(PEPTIDE, alleles=",".join(LARGE_POOL), cls="mhc1")
    print(f"binder_score over {len(LARGE_POOL)} alleles: {time.time() - _t0:.1f} s")

    _rows = "\n".join(
        f"| {i} | {esc(b.allele)} | {b.binder_rank:.2f} | {b.band} | {b.p_binder:.4f} |"
        for i, b in enumerate(large, 1)
    )
    mo.md(
        f"""
    **`{PEPTIDE}` against a {len(LARGE_POOL)}-allele pool**

    | position | allele | binder_rank | band | p_binder |
    |--:|---|--:|---|--:|
    {_rows}
    """
    )
    return (large,)


@app.cell
def _(SMALL_POOL, esc, large, mo, small):
    _big = {b.allele: b for b in large}
    _small_pos = {b.allele: i for i, b in enumerate(small, 1)}
    _large_pos = {b.allele: i for i, b in enumerate(large, 1)}

    _rows = "\n".join(
        f"| {esc(a)} | {_small_pos[a]} | {_large_pos[a]} | "
        f"{next(b for b in small if b.allele == a).binder_rank:.2f} | {_big[a].binder_rank:.2f} | "
        f"{next(b for b in small if b.allele == a).p_binder:.4f} | {_big[a].p_binder:.4f} |"
        for a in SMALL_POOL
    )
    mo.md(
        f"""
    **The shared alleles, side by side**

    | allele | position in 3-pool | position in {len(large)}-pool | binder_rank in 3-pool | binder_rank in {len(large)}-pool | p_binder in 3-pool | p_binder in {len(large)}-pool |
    |---|--:|--:|--:|--:|--:|--:|
    {_rows}

    HLA-A\\*02:01 moves from position {_small_pos["HLA-A*02:01"]} to position
    {_large_pos["HLA-A*02:01"]} — three A\\*02 relatives out-score it once they are in the pool — while
    its `binder_rank` and `p_binder` are byte-identical in both runs. That is the design: both numbers
    are calibrated against a fixed random-peptide background *per allele*, so neither depends on which
    other alleles you happened to ask about.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1.5 Where the pool *does* bite: ranking candidates

    The pool-dependence that actually costs you is one level up. A pipeline that says "take the top 15%
    of candidates" is computing a percentile **within the candidate list**, and that renormalises every
    time the list changes.

    Below, the same peptide `NLVPMVATV` sits in two 7-peptide candidate lists scored on the same allele
    (`HLA-A*02:01`):

    * a **mixed** list — one A\*02:01 epitope among six epitopes restricted to other alleles;
    * an **A\*02:01** list — seven epitopes all restricted to A\*02:01.

    Watch the within-list percentile against `p_binder`.
    """)
    return


@app.cell
def _(mo, store, time):
    ALLELE = "HLA-A*02:01"
    MIXED = ["NLVPMVATV", "TPRVTGGGAM", "RAKFKQLL", "AVFDRKSDAK", "LVVDFSQFSR", "KRWIILGLNK", "IPSINVHHY"]
    A2_ONLY = ["NLVPMVATV", "GILGFVFTL", "GLCTLVAML", "YVLDHLIVV", "LLWNGPMAV", "LLLDRLNQL", "ELAGIGILTV"]


    def score_pool(peptides):
        """binder_rank / p_binder for each peptide on ALLELE, best-first."""
        out = [store.binder_score(p, alleles=ALLELE, cls="mhc1")[0] for p in peptides]
        out.sort(key=lambda b: b.binder_rank)
        return out


    _t0 = time.time()
    mixed_scored = score_pool(MIXED)
    a2_scored = score_pool(A2_ONLY)
    print(f"14 peptide-allele scores on a warm calibrator: {time.time() - _t0:.1f} s")


    def _table(scored):
        return "\n".join(
            f"| {i} | {b.peptide} | {100.0 * i / len(scored):.0f}% | {b.binder_rank:.2f} | {b.p_binder:.4f} |"
            for i, b in enumerate(scored, 1)
        )


    mo.md(
        f"""
    **Mixed candidate list**

    | position | peptide | within-list percentile | binder_rank | p_binder |
    |--:|---|--:|--:|--:|
    {_table(mixed_scored)}

    **A\\*02:01-only candidate list**

    | position | peptide | within-list percentile | binder_rank | p_binder |
    |--:|---|--:|--:|--:|
    {_table(a2_scored)}
    """
    )
    return a2_scored, mixed_scored


@app.cell
def _(a2_scored, mixed_scored, mo):
    _m = [b.peptide for b in mixed_scored].index("NLVPMVATV") + 1
    _a = [b.peptide for b in a2_scored].index("NLVPMVATV") + 1
    _b = next(b for b in mixed_scored if b.peptide == "NLVPMVATV")

    mo.md(
        f"""
    ### The point

    | quantity for `NLVPMVATV` on HLA-A\\*02:01 | mixed list | A\\*02:01-only list |
    |---|--:|--:|
    | position in list | **{_m}** | {_a} |
    | within-list percentile | **{100.0 * _m / len(mixed_scored):.0f}%** | {100.0 * _a / len(a2_scored):.0f}% |
    | `binder_rank` | **{_b.binder_rank:.2f}** | **{_b.binder_rank:.2f}** |
    | `p_binder` | **{_b.p_binder:.4f}** | **{_b.p_binder:.4f}** |

    Nothing about the peptide changed. Only the company it keeps changed, and the within-list
    percentile moved by {abs(100.0 * _a / len(a2_scored) - 100.0 * _m / len(mixed_scored)):.0f}
    percentage points. `binder_rank` and `p_binder` did not move at all.

    **Use `binder_rank` to sort a list. Use `p_binder` whenever the number has to mean the same thing
    outside the list it was computed in** — an absolute cut-off, a feature fed to a downstream model, a
    comparison against a run made last month with a different candidate set.
    """
    )
    return


@app.cell
def _(mhcmatch, mo):
    mo.md(
        f"""
    ---

    ### Also available from the same `Store`

    ```python
    store.is_binder("NLVPMVATV", "HLA-A*02:01")
    store.scan_protein(my_protein_seq, cls="mhc1", correction="bh")   # FDR-controlled windows
    store.anchor_model("mhc1", background="proteome")                 # the presentation scorer
    store.affinity_model("mhc1").amplitude("NLVPMVATL", "NLVPMVATV", "HLA-A*02:01")
    ```

    and on the command line:

    ```bash
    mhcmatch binder NLVPMVATV --alleles 'HLA-A*02:01,HLA-B*07:02' --cls mhc1
    ```

    mhcmatch {mhcmatch.__version__ if hasattr(mhcmatch, "__version__") else ""}
    """
    )
    return


if __name__ == "__main__":
    app.run()
