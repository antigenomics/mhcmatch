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
    # 2 — Physicochemical epitope features (`mhcmatch.immuno`)

    **What this demonstrates.** `mhcmatch.immuno.features()` turns a peptide into **141 numbers**:
    `length`, plus 20 amino-acid scales x 7 statistics. The two design choices that make it more than a
    bag of descriptors are both parameters you can see and change here:

    1. **Which positions count as TCR-facing** — `ANCHOR_SCHEMES` keeps every anchor definition in the
       toolchain selectable, plus a continuous `"contact"` weighting derived from crystal contacts.
    2. **How residues are aggregated** — sums and means are the established descriptors, but a
       *contiguous* hydrophobic stretch is a different object from the same residues scattered, and no
       sum can express that. Hence `run_max` / `run_n` / `run_frac`.

    **What you should conclude.** The run statistics carry information that is provably absent from the
    sums: section 2.3 shows three peptides with **identical amino-acid composition** and therefore
    identical `*_sum`, `*_mean`, `*_min`, `*_max` and `*_run_frac`, but different `run_max` / `run_n`.
    And section 2.4 shows that the contact profile recovers an anchor set *unsupervised* which matches
    neither shipped scheme.

    Nothing in this notebook needs the reference panel — the scales and the contact profile are
    vendored in the package, so this runs offline in under a second.
    """)
    return


@app.cell
def _():
    from statistics import median

    from mhcmatch import immuno
    from mhcmatch.data import aa_tables

    EPITOPE = "GILGFVFTL"  # influenza A M1 58-66, HLA-A*02:01

    _s = immuno.scales()
    print(f"scale families requested: {immuno.DEFAULT_SCALES}")
    print(f"-> {len(_s)} scales: {sorted(_s)}")
    print(f"anchor schemes: {sorted(immuno.ANCHOR_SCHEMES)} + 'contact'")
    print(f"features({EPITOPE!r}) -> {len(immuno.features(EPITOPE))} values")
    print(f"  = 1 length + {len(_s)} scales x 7 statistics")
    return EPITOPE, aa_tables, immuno, median


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.1 The three anchor schemes, plus the contact scheme

    `position_weights` returns one weight per position: `0` at masked (anchor) positions, `1` at
    TCR-facing ones. `"contact"` returns a *continuous* weight instead, rescaled so the surviving
    positions have mean 1 — which puts it on the same footing as the binary schemes.

    The three binary schemes are the three incompatible class-I anchor definitions that coexist in this
    toolchain; they are kept selectable rather than collapsed into a constant.
    """)
    return


@app.cell
def _(EPITOPE, immuno, mo):
    _cp = immuno.contact_profile("mhc1")
    _schemes = {
        "full": immuno.position_weights(EPITOPE, "mhc1", "full"),
        "p2_pomega": immuno.position_weights(EPITOPE, "mhc1", "p2_pomega"),
        "pockets": immuno.position_weights(EPITOPE, "mhc1", "pockets"),
        "contact": immuno.position_weights(EPITOPE, "mhc1", "contact", contact_profile=_cp),
    }

    _hdr = " | ".join(f"P{i + 1} ({c})" for i, c in enumerate(EPITOPE))
    _sep = "|".join(["--:"] * len(EPITOPE))
    _rows = "\n".join(
        f"| `{name}` | " + " | ".join(f"{w:.2f}" for w in ws) + f" | {sum(1 for w in ws if w == 0)} |"
        for name, ws in _schemes.items()
    )
    mo.md(
        f"""
    **Per-position weights for `{EPITOPE}`**

    | scheme | {_hdr} | n masked |
    |---|{_sep}|--:|
    {_rows}
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.2 What a masked position does to the features

    Masking is not a cosmetic relabelling — it changes the weighted sums, and it changes the run
    structure, which is the point of section 2.3.
    """)
    return


@app.cell
def _(EPITOPE, immuno, mo):
    _cp2 = immuno.contact_profile("mhc1")
    _cols = ["MJ_sum", "MJ_mean", "KyteDoolittle_sum", "KyteDoolittle_mean", "KF1_sum", "KF6_sum"]
    _variants = {
        "full": immuno.features(EPITOPE, scheme="full"),
        "p2_pomega": immuno.features(EPITOPE, scheme="p2_pomega"),
        "pockets": immuno.features(EPITOPE, scheme="pockets"),
        "contact": immuno.features(EPITOPE, scheme="contact", contact_profile=_cp2),
    }
    _rows = "\n".join(
        f"| `{name}` | " + " | ".join(f"{f[c]:.3f}" for c in _cols) + " |" for name, f in _variants.items()
    )
    mo.md(
        f"""
    **Six of the 141 features for `{EPITOPE}`, under each scheme**

    | scheme | {" | ".join(f"`{c}`" for c in _cols)} |
    |---|{"|".join(["--:"] * len(_cols))}|
    {_rows}
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.3 Run statistics express contiguity — and no sum can

    Take the CMV epitope `NLVPMVATV` and two **permutations** of it. All three have exactly the same
    amino-acid composition, so every composition-derived statistic is identical by construction:
    `sum`, `mean`, `min`, `max`, and even `run_frac` (which counts *how many* positions are above the
    scale's median, not where they are).

    * `clustered` puts every above-median Kyte-Doolittle residue together;
    * `interleaved` alternates them with the below-median ones;
    * `native` is the real epitope.

    Only `run_max` and `run_n` can tell them apart.
    """)
    return


@app.cell
def _(aa_tables, immuno, median, mo):
    _kd = aa_tables.HYDROPHOBICITY["KyteDoolittle"]
    _thr = median(_kd.values())

    NATIVE = "NLVPMVATV"
    _hi = sorted(c for c in NATIVE if _kd[c] > _thr)
    _lo = sorted(c for c in NATIVE if _kd[c] <= _thr)
    CLUSTERED = "".join(_hi + _lo)
    _h, _l, _mix = list(_hi), list(_lo), []
    while _h or _l:
        if _h:
            _mix.append(_h.pop(0))
        if _l:
            _mix.append(_l.pop(0))
    INTERLEAVED = "".join(_mix)

    assert sorted(NATIVE) == sorted(CLUSTERED) == sorted(INTERLEAVED)
    print(f"Kyte-Doolittle median threshold: {_thr}")
    print(f"native      {NATIVE}")
    print(f"clustered   {CLUSTERED}")
    print(f"interleaved {INTERLEAVED}")
    print("identical composition:", sorted(NATIVE) == sorted(INTERLEAVED))

    _perms = {"native": NATIVE, "clustered": CLUSTERED, "interleaved": INTERLEAVED}
    _f = {k: immuno.features(v, scheme="full") for k, v in _perms.items()}
    _rows = "\n".join(
        f"| {k} | `{v}` | {_f[k]['KyteDoolittle_sum']:.2f} | {_f[k]['KyteDoolittle_mean']:.3f} | "
        f"{_f[k]['KyteDoolittle_run_frac']:.3f} | {_f[k]['KyteDoolittle_run_max']:.0f} | "
        f"{_f[k]['KyteDoolittle_run_n']:.0f} |"
        for k, v in _perms.items()
    )
    mo.md(
        f"""
    **Kyte-Doolittle statistics, `scheme="full"` (no masking, so the sums are pure composition)**

    | arrangement | peptide | `_sum` | `_mean` | `_run_frac` | `_run_max` | `_run_n` |
    |---|---|--:|--:|--:|--:|--:|
    {_rows}

    The first three columns are constant across the row set — they *cannot* separate these peptides.
    `run_max` and `run_n` do, and they say something structural: `clustered` presents one unbroken
    7-residue hydrophobic face, `interleaved` presents three short ones.
    """
    )
    return (CLUSTERED,)


@app.cell
def _(CLUSTERED, immuno, mo):
    _rows = "\n".join(
        f"| `{s}` | {f['KyteDoolittle_sum']:.2f} | {f['KyteDoolittle_run_max']:.0f} | {f['KyteDoolittle_run_n']:.0f} |"
        for s, f in (
            (s, immuno.features(CLUSTERED, scheme=s)) for s in ("full", "p2_pomega", "pockets")
        )
    )
    mo.md(
        f"""
    ### A masked anchor **breaks** a run rather than bridging it

    Same peptide (`{CLUSTERED}`), different masking. A buried residue sitting between two exposed
    hydrophobics does not make them contiguous from the TCR's point of view, so the run walker treats a
    zero-weight position as a break — not as an absence.

    | scheme | `KyteDoolittle_sum` | `_run_max` | `_run_n` |
    |---|--:|--:|--:|
    {_rows}

    One unbroken run of 7 becomes two runs, the longest of 5, as soon as P2 is masked.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.4 The contact profile

    `immuno.contact_profile("mhc1")` is backed by 8,062 TCR-peptide residue contacts over 370 crystal
    structures. It turns a contact *frequency* into a *weight* in two derived steps, neither tuned:

    * **zeroing** — a position contacted less than half as often as a uniform footprint would predict
      (`frac < 1/(2L)`) is treated as not TCR-facing;
    * **rescaling** — the survivors are scaled to mean 1, so the weighted statistics land on the same
      scale as the binary schemes.

    On class-I 9-mers this zeroes exactly **P1, P2, P3 and P-omega** — an anchor set recovered from the
    contact data alone, without being told that anchors exist.
    """)
    return


@app.cell
def _(immuno, mo):
    _p1 = immuno.contact_profile("mhc1")
    _w9 = _p1(9)
    _zeroed = [i + 1 for i, w in enumerate(_w9) if w == 0.0]
    _nz = [w for w in _w9 if w > 0]

    _rows = "\n".join(
        f"| P{i + 1} | {w:.3f} | {'masked' if w == 0.0 else 'TCR-facing'} |" for i, w in enumerate(_w9)
    )
    print(f"zeroed positions (1-based): {_zeroed}")
    print(f"mean weight of surviving positions: {sum(_nz) / len(_nz):.6f}")
    print(f"class-I 12-mer (a length with thin structural support): {[round(w, 2) for w in _p1(12)]}")
    print(f"class-II 15-mer: {[round(w, 2) for w in immuno.contact_profile('mhc2')(15)]}")

    mo.md(
        f"""
    **`contact_profile("mhc1")(9)`**

    | position | weight | call |
    |---|--:|---|
    {_rows}

    Zeroed: **P{", P".join(str(z) for z in _zeroed[:-1])} and P-omega (P{_zeroed[-1]})**.

    That set is **neither** shipped scheme: `p2_pomega` misses P1 and P3, and `pockets` additionally
    masks P8 — which the contact data ranks among the most-contacted positions. Which definition wins
    is therefore an ablation with a reported number, not a constant to hardcode.

    Lengths with too few observed structures fall back to the class's pooled *relative*-position
    profile, interpolated onto the requested length: the TCR footprint scales with the peptide, so
    relative position transfers where absolute position does not.
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.5 Building a feature matrix

    `feature_names()` gives the column order without a dict round-trip, so a matrix is one
    comprehension. `length` is emitted deliberately: ligand length distribution is allele-specific and
    is part of what distinguishes a real ligand set — it is signal here, not a nuisance to regress out.
    """)
    return


@app.cell
def _(immuno):
    PEPTIDES = ["GILGFVFTL", "NLVPMVATV", "GLCTLVAML", "TPRVTGGGAM", "RAKFKQLL", "AVFDRKSDAK"]

    _cp3 = immuno.contact_profile("mhc1")
    cols = immuno.feature_names()
    matrix = [
        [immuno.features(p, scheme="contact", contact_profile=_cp3)[c] for c in cols] for p in PEPTIDES
    ]

    print(f"matrix: {len(matrix)} peptides x {len(cols)} features")
    print(f"first 8 columns: {cols[:8]}")
    for _p, _row in zip(PEPTIDES, matrix):
        print(f"  {_p:<11s} length={_row[0]:.0f}  KF1_sum={_row[1]:+.3f}  KF1_run_max={_row[5]:.0f}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---

    ### Notes

    * **Class II ignores `scheme`.** The class-II core P1/P4/P6/P9 definition *is* agreed across the
      toolchain, so `position_weights(..., cls="mhc2")` always masks the register-anchored core. Pass
      `register_start=` from `AnchorModel.best_register` so the annotated frame matches the scored one;
      `None` uses the allele-agnostic heuristic register.
    * **Non-standard residues are dropped, not zeroed.** A zero is a real value on a centred scale like
      Kidera, so scoring `X` as 0 would be a silent bias.
    * **Any AAindex-style table works.** All scales are plain `dict[str, float]`, so
      `features(..., scale_names=(...))` accepts anything in `mhcmatch.data.aa_tables`.
    * `python -m mhcmatch.immuno` runs the module's own self-check against published constants.
    """)
    return


if __name__ == "__main__":
    app.run()
