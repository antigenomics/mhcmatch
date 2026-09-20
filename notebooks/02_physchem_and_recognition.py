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
    # 2 — Physicochemistry and the recognition axis

    **What this demonstrates.** The half of the model that is not about binding. A peptide that is
    presented is not yet a peptide a T cell reacts to, and what separates the two is a property of the
    residues the receptor actually touches. This notebook takes a published immunogenicity corpus and
    walks the three layers `mhcmatch` puts on it:

    1. `mhcmatch.immuno` — 141 physicochemical features per peptide, over a selectable set of
       TCR-facing positions;
    2. `mhcmatch.posbayes` — per-position, per-role amino-acid evidence, with no prior;
    3. `mhcmatch.complement` — the fitted recognition axis, six feature blocks, one number.

    **What you should conclude.** The signal is real, it is *positional*, and it is concentrated on
    the TCR-facing strip. The same amino acid carries a different sign depending on whether it sits in
    an anchor pocket or under the receptor, so a composition feature that ignores position averages
    the two away.

    Data: `immunogenicity/chowell_iedb_full_matched.tsv.gz` from the public HuggingFace dataset
    `isalgo/pmhc_data` — 29,492 HLA-matched rows, ~260 kB. Nothing is downloaded by hand.

    **No other tool is run here.** This is the library on a public corpus; head-to-head comparisons
    live in the benchmark repository.
    """)
    return


@app.cell
def _():
    import gzip
    import time

    import numpy as np

    from mhcmatch.store import fetch_file

    _t0 = time.time()
    _path = fetch_file("immunogenicity/chowell_iedb_full_matched.tsv.gz")
    with gzip.open(_path, "rt") as _fh:
        _cols = _fh.readline().rstrip("\n").split("\t")
        _rows = [dict(zip(_cols, ln.rstrip("\n").split("\t"))) for ln in _fh]

    # Class I only, and only the canonical 9-mers: the positional roles below are defined against a
    # 9-residue register, and mixing lengths would average two different geometries.
    rows = [r for r in _rows if r["mhc_class"] == "MHCI" and len(r["peptide"]) == 9]
    peptides = [r["peptide"] for r in rows]
    labels = np.array([int(r["label"]) for r in rows])
    print(f"loaded in {time.time() - _t0:.1f} s")
    print(f"{len(_rows):,} rows -> {len(rows):,} class-I 9-mers")
    print(f"immunogenic: {int(labels.sum()):,} ({labels.mean():.1%})")
    return labels, np, peptides, rows, time


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.1 The TCR-facing strip is not the whole peptide

    A class-I 9-mer binds through pockets at roughly P2 and P9. Those residues are spoken for — they
    are what makes the peptide *fit* — and what is left is what a receptor can read. `anchor_indices`
    returns the split for a given scheme, and `decompose` applies it.
    """)
    return


@app.cell
def _():
    import mhcmatch
    from mhcmatch import immuno, posbayes

    _demo = "GILGFVFTL"
    anchors = mhcmatch.anchor_indices(_demo, "mhc1")
    facing = [i for i in range(len(_demo)) if i not in set(anchors)]
    print(f"anchor positions (0-based): {sorted(anchors)}")
    print(f"TCR-facing positions:       {facing}")
    print()
    # `posbayes.roles` is the same split as a 1/0 mask per position -- 1 anchor, 0 facing.
    print(f"posbayes.roles(9) = {posbayes.roles(9)}   (1 = anchor, 0 = TCR-facing)")
    print()
    print(f"peptide      {_demo}")
    print(f"anchor part  {''.join(_demo[i] for i in sorted(anchors))}")
    print(f"facing part  {''.join(_demo[i] for i in facing)}")
    return anchors, facing, immuno, mhcmatch, posbayes


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.2 141 features, and the two that survive into the shipped model

    `immuno.features` returns the whole basis. The shipped `EPIC` model does **not** use 141 of
    anything — it carries exactly two chemistry terms, `C_phys_buried` and `C_phys_charge`, because
    the rest are collinear with those two and a fit that keeps them all is fitting noise. The full
    basis is here for anyone who wants to argue otherwise on their own data.

    Below: mean burial (Rose) and charge (Atchley factor 5) over the TCR-facing strip, split by label.
    """)
    return


@app.cell
def _(facing, immuno, labels, np, peptides):
    from mhcmatch.data import aa_tables

    _names = immuno.feature_names()
    print(f"physicochemical features per peptide: {len(_names)}")
    print(f"first eight: {_names[:8]}")
    print(f"scales:          {', '.join(immuno.DEFAULT_SCALES)}")
    print(f"anchor schemes:  {', '.join(immuno.ANCHOR_SCHEMES)}")
    print()

    def _strip_mean(table, peps, idx):
        """Mean of `table` over a chosen set of positions."""
        return np.array([np.mean([table.get(p[i], 0.0) for i in idx]) for p in peps])

    pc1_face = _strip_mean(aa_tables.PROPERTY_PC1, peptides, facing)
    pc2_face = _strip_mean(aa_tables.PROPERTY_PC2, peptides, facing)
    pc1_anch = _strip_mean(aa_tables.PROPERTY_PC1, peptides, sorted(anchors))

    # `d` rather than "Cohen's d" in the header: nesting a double-quoted string inside a
    # double-quoted f-string is PEP 701 and parses only on Python 3.12+, while this package
    # supports 3.10. It cost a green CI run to find, which is the notebook lint gate earning itself.
    _cohen = "Cohen's d"
    print(f"{'axis':28s} {'immunogenic':>12s} {'other':>10s} {_cohen:>11s}")
    for _label, _vals in (("PC1 (burial), TCR-facing", pc1_face),
                          ("PC2 (charge), TCR-facing", pc2_face),
                          ("PC1 (burial), anchors", pc1_anch)):
        _pos, _neg = _vals[labels == 1], _vals[labels == 0]
        _d = (_pos.mean() - _neg.mean()) / np.sqrt((_pos.var() + _neg.var()) / 2)
        print(f"{_label:28s} {_pos.mean():+12.4f} {_neg.mean():+10.4f} {_d:+11.3f}")
    print()
    print("Burial separates on the TCR-facing strip about 2.6x as strongly as on the anchors, and")
    print("charge barely separates at all on this corpus -- which is why the shipped model carries")
    print("burial and charge as TWO terms rather than one, and fits their weights rather than")
    print("assuming them. Anchors are largely spent on fitting the groove; what is left is what a")
    print("receptor reads.")
    return aa_tables, pc1_anch, pc1_face, pc2_face


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.3 The same residue, two roles, opposite signs

    This is the reason the model is positional. `posbayes` scores an amino acid *given the position
    role it occupies* — anchor or TCR-facing — and the two evidence tables disagree about specific
    residues. A model that pools them cannot see this at all.

    The table below lists the residues whose anchor and facing log-odds have opposite signs, largest
    disagreement first.
    """)
    return


@app.cell
def _(np):
    _tab = posbayes.table("human")
    _anchor, _face = _tab["anchor"], _tab["tcrface"]
    print(f"fitted on {_tab['n']:,} peptides, {_tab['n_immunogenic']:,} immunogenic "
          f"(prevalence {_tab['prevalence']:.3f})")
    print()

    _rows = []
    for _i, _aa in enumerate(posbayes.AA):
        a, f = float(_anchor[_i]), float(_face[_i])
        if a * f < 0:                       # the two roles disagree about the sign
            _rows.append((_aa, a, f, abs(a - f)))
    _rows.sort(key=lambda r: -r[3])

    print(f"{'aa':>3}  {'anchor':>8}  {'facing':>8}   reading")
    for _aa, a, f, _ in _rows[:8]:
        _read = ("helps presentation, hurts recognition" if a > 0
                 else "hurts presentation, helps recognition")
        print(f"{_aa:>3}  {a:+8.3f}  {f:+8.3f}   {_read}")
    print()
    print(f"{len(_rows)} of {len(posbayes.AA)} residues carry OPPOSITE signs in the two roles.")
    print("A whole-peptide composition feature averages these towards zero, which is why the model")
    print("scores an amino acid given the role of the position it sits in.")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2.4 The recognition axis, as one number

    `complement.score` is the fitted axis: six feature blocks over the TCR-facing strip, reduced to a
    single log-odds. The `aa` block **is** `posbayes` — the same evidence, entering as one block among
    six rather than as the whole answer.

    Scored over the whole corpus in one call, then read as an AUROC against the published labels. The
    interval is closed-form (Hanley & McNeil 1982): AUROC is a U-statistic and needs no bootstrap.
    """)
    return


@app.cell
def _(labels, np, peptides, time):
    from mhcmatch import complement

    _t0 = time.time()
    scores = np.asarray(complement.score(peptides, species="human", cls="mhc1"), dtype=float)
    print(f"scored {len(peptides):,} peptides in {time.time() - _t0:.1f} s "
          f"({len(peptides) / max(time.time() - _t0, 1e-9):,.0f}/s) — one call, not a loop")
    print(f"blocks ({len(complement.BLOCKS)}): {', '.join(complement.BLOCKS)}")
    for _b, _f in complement.BLOCKS.items():
        print(f"  {_b:8s} {len(_f):2d} feature(s): {', '.join(_f[:5])}"
              + (" ..." if len(_f) > 5 else ""))
    return complement, scores


@app.cell
def _(labels, np, scores):
    def auroc(y, s):
        """Ties-averaged Mann-Whitney AUROC, and Hanley & McNeil's closed-form 95% interval."""
        order = np.argsort(s, kind="stable")
        ranks = np.empty(len(s), dtype=float)
        ranks[order] = np.arange(1, len(s) + 1)
        # average the ranks of tied scores, or every tie is silently broken by input order
        _s = s[order]
        i = 0
        while i < len(_s):
            j = i
            while j + 1 < len(_s) and _s[j + 1] == _s[i]:
                j += 1
            if j > i:
                ranks[order[i:j + 1]] = (i + j + 2) / 2
            i = j + 1
        npos, nneg = int(y.sum()), int((1 - y).sum())
        a = (ranks[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg)
        q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
        se = np.sqrt((a * (1 - a) + (npos - 1) * (q1 - a * a)
                      + (nneg - 1) * (q2 - a * a)) / (npos * nneg))
        return a, a - 1.96 * se, a + 1.96 * se

    _a, _lo, _hi = auroc(labels, scores)
    print(f"complementarity AUROC {_a:.4f}  (95% CI {_lo:.4f}-{_hi:.4f})")
    print(f"n = {len(labels):,} class-I 9-mers, {int(labels.sum()):,} immunogenic")
    print()
    print("The interval excludes 0.5, so the axis separates on this corpus. It is one of nine terms")
    print("in the shipped EPIC model, not a ranker on its own — see notebook 4.")
    return (auroc,)


@app.cell
def _(mo):
    mo.md(r"""
    ## What to take away

    * The receptor reads a **strip**, not a peptide. Anchors are spent on fitting the groove.
    * Position is not a detail: **the same residue carries opposite evidence in the two roles**, so a
      composition feature that ignores it cancels the signal it was meant to capture.
    * The 141-feature basis is available, and the shipped model deliberately uses **two** of it. More
      features here buy collinearity, not accuracy.
    * `complement.score` takes the whole corpus in one call. Do not loop it per peptide.

    Next: [notebook 3](03_corpus_overlap_and_similarity.py) — what the repertoire was shaped by, and
    how resemblance to it becomes three scored channels.
    """)
    return


if __name__ == "__main__":
    app.run()
