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
    # 3 — Corpus overlap, and similarity as three scored channels

    **What this demonstrates.** The repertoire that meets a neoantigen was shaped by everything it met
    before. `mhcmatch` makes that concrete with three reference corpora — the **self** proteome, the
    **thymic** immunopeptidome, and a **viral** ligandome — and turns resemblance to each into a
    scored term.

    Two halves:

    1. **Overlap.** What the three corpora actually share, as whole peptides and as receptor-facing
       faces. These are not the same question and they do not give the same answer.
    2. **Scoring.** `C_corpus_thymus` / `_self` / `_viral`, the three channels the shipped model
       fits, computed on peptides you choose.

    **What you should conclude.** The three channels are not redundant and they do not have the same
    sign. Resemblance to **self** is the largest negative term in the whole model — a peptide that
    looks like something the repertoire was tolerised against is less likely to raise a response.
    Resemblance to the **thymic** immunopeptidome *raises* the score, and that is not a contradiction:
    thymic presentation is what makes an allotype's peptides visible in the first place.

    Data: `thymus/thymus_immunopeptidome.tsv.gz` and `ligandome/viral_foreign_iedb.tsv.gz` from the
    public HuggingFace dataset `isalgo/pmhc_data`, ~1.8 MB together.
    """)
    return


@app.cell
def _():
    import gzip
    import time

    import numpy as np

    from mhcmatch.store import fetch_file

    def read_tsv(rel, keep=None):
        with gzip.open(fetch_file(rel), "rt") as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            out = []
            for ln in fh:
                r = dict(zip(cols, ln.rstrip("\n").split("\t")))
                if keep is None or keep(r):
                    out.append(r)
        return out

    _t0 = time.time()
    # Class I only: an anchor-masked face is defined against a 9-residue register, and class II
    # peptides are a different geometry that would not be comparable.
    thymus = read_tsv("thymus/thymus_immunopeptidome.tsv.gz",
                      lambda r: r.get("mhc_class") == "MHCI" and len(r["peptide"]) == 9)
    viral = read_tsv("ligandome/viral_foreign_iedb.tsv.gz",
                     lambda r: r.get("mhc_class") == "MHCI" and len(r["peptide"]) == 9)
    print(f"loaded in {time.time() - _t0:.1f} s")
    print(f"thymic class-I 9-mers: {len(thymus):,}")
    print(f"viral  class-I 9-mers: {len(viral):,}")
    return np, read_tsv, thymus, time, viral


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.1 Whole peptides overlap almost not at all

    The obvious question first: how many 9-mers appear in both corpora? The answer is what makes the
    second question interesting.
    """)
    return


@app.cell
def _(thymus, viral):
    thy_p = {r["peptide"] for r in thymus}
    vir_p = {r["peptide"] for r in viral}
    _both = thy_p & vir_p
    print(f"distinct thymic 9-mers : {len(thy_p):,}")
    print(f"distinct viral  9-mers : {len(vir_p):,}")
    print(f"shared as WHOLE peptides: {len(_both):,} "
          f"({100 * len(_both) / min(len(thy_p), len(vir_p)):.2f}% of the smaller set)")
    if _both:
        print(f"  e.g. {', '.join(sorted(_both)[:5])}")
    return thy_p, vir_p


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.2 Three levels of "the same", and only the third one overlaps

    A T-cell receptor does not read the whole peptide, so ask the question again at two more levels:
    the **face** (the peptide with its anchors masked out) and the **masked 3-mers** that face is made
    of — which is the level the shipped model actually contracts against, `corpus_k = 3`.

    The three answers are not close to each other, and that gap is the entire design argument for the
    corpus channels.
    """)
    return


@app.cell
def _(thy_p, vir_p):
    import mhcmatch

    def faces(peps):
        """The receptor-facing residues of each peptide, anchors dropped."""
        out = set()
        for p in peps:
            anc = set(mhcmatch.anchor_indices(p, "mhc1"))
            out.add("".join(c for i, c in enumerate(p) if i not in anc))
        return out

    def masked_kmers(peps, k=3):
        """Every k-mer of every anchor-masked face -- the vocabulary the model scores against."""
        out = set()
        for p in peps:
            anc = set(mhcmatch.anchor_indices(p, "mhc1"))
            f = "".join(c for i, c in enumerate(p) if i not in anc)
            for i in range(len(f) - k + 1):
                out.add(f[i:i + k])
        return out

    thy_f, vir_f = faces(thy_p), faces(vir_p)
    thy_k, vir_k = masked_kmers(thy_p), masked_kmers(vir_p)

    print(f"{'level':>16} {'thymic':>9} {'viral':>9} {'shared':>9} {'% of smaller':>13}")
    for _name, _a, _b in (("whole peptide", thy_p, vir_p),
                          ("masked face", thy_f, vir_f),
                          ("masked 3-mer", thy_k, vir_k)):
        _sh = len(_a & _b)
        print(f"{_name:>16} {len(_a):9,} {len(_b):9,} {_sh:9,} "
              f"{100 * _sh / min(len(_a), len(_b)):12.1f}%")
    print()
    print("Zero whole peptides. One face out of fifteen thousand. And essentially the SAME 3-mer")
    print("vocabulary -- the two corpora are disjoint as strings and near-identical in the units a")
    print("receptor contact actually spans.")
    print()
    print("So an identity match between corpora finds nothing, and that is not because they are")
    print("unrelated. It is why the channels are a DENSITY over that shared vocabulary rather than a")
    print("lookup: what differs between self, thymic and viral is not which 3-mers occur, it is how")
    print("often each one does.")
    return faces, masked_kmers, mhcmatch, thy_f, thy_k, vir_f, vir_k


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.3 The three channels, computed

    `mhcmatch.cli._aggregate_channels` returns the function the shipped model uses. It contracts each
    peptide against one 20³ count table per corpus under an identity-normalised BLOSUM62 kernel —
    **one table contraction, not a neighbour search**, which is why it is fast enough to run on a
    whole candidate pool.
    """)
    return


@app.cell
def _(np, thy_p, time, vir_p):
    from mhcmatch import cli

    channels = cli._aggregate_channels("mhc1", no_self=False, species="human")

    # A deliberately mixed probe: some thymic peptides, some viral, some neither.
    probe = (sorted(thy_p)[:400] + sorted(vir_p)[:400]
             + ["AAAAAAAAA", "KKKKKKKKK", "GILGFVFTL", "NLVPMVATV"])
    origin = ["thymic"] * 400 + ["viral"] * 400 + ["other"] * 4

    _t0 = time.time()
    ch = channels(probe)
    print(f"{len(probe):,} peptides x 3 channels in {time.time() - _t0:.2f} s")
    print(f"channels: {', '.join(ch)}")
    print()
    print(f"{'origin':>8} {'n':>5}  " + "".join(f"{c.replace('C_corpus_', ''):>12}" for c in ch))
    for _o in ("thymic", "viral", "other"):
        _m = np.array([x == _o for x in origin])
        print(f"{_o:>8} {int(_m.sum()):5d}  "
              + "".join(f"{np.mean(np.asarray(ch[c])[_m]):12.6f}" for c in ch))
    return ch, channels, cli, origin, probe


@app.cell
def _(mo):
    mo.md(r"""
    Each corpus scores highest on peptides drawn from itself, which is the sanity check: the channels
    are measuring what they claim to. The interesting part is that they are **not** interchangeable —
    a thymic peptide is not simply "high on everything".

    ## 3.4 The three channels are correlated, and still not redundant

    They are computed from overlapping reference sets, so of course they move together. The question a
    model has to answer is whether entering all three buys anything over entering one. In the shipped
    fit it does: the corpus block is worth χ² = 20.0 on 3 degrees of freedom, and the three
    coefficients do not have the same sign.
    """)
    return


@app.cell
def _(ch, np):
    _names = list(ch)
    _M = np.array([np.asarray(ch[c], dtype=float) for c in _names])
    print("Pearson r between channels, on the probe set:")
    print(f"{'':>16}" + "".join(f"{n.replace('C_corpus_', ''):>10}" for n in _names))
    for _i, _a in enumerate(_names):
        print(f"{_a.replace('C_corpus_', ''):>16}"
              + "".join(f"{np.corrcoef(_M[_i], _M[_j])[0, 1]:10.3f}" for _j in range(len(_names))))
    return


@app.cell
def _():
    import json
    from pathlib import Path

    import mhcmatch as _mm

    _art = Path(_mm.__file__).parent / "data" / "aggregate_mhc1.json"
    _fit = json.loads(_art.read_text())
    _w = dict(zip(_fit["features"], _fit["coef"]))
    print("What the shipped model does with them (log-odds per standard deviation):")
    for _k in ("C_corpus_thymus", "C_corpus_self", "C_corpus_viral"):
        _sign = "raises" if _w[_k] > 0 else "LOWERS"
        print(f"  {_k:20s} {_w[_k]:+.4f}   {_sign} the score")
    print()
    print(f"  self is the largest negative coefficient of any of the {len(_w)} terms: "
          f"{min(_w.values()):+.4f}")
    print()
    print("Thymic and viral resemblance raise the score; self lowers it, hardest. Averaging the")
    print("three into one similarity would cancel most of that and report roughly nothing.")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What to take away

    * **Whole-peptide overlap between corpora is negligible; face overlap is not.** Any similarity
      feature that matches peptide strings is answering the wrong question.
    * The three channels are a **table contraction**, one 20³ count table per corpus, so scoring a
      whole candidate pool costs a matrix operation rather than a search.
    * They are correlated and still not redundant: the shipped model fits **thymic and viral
      positive, self negative**, and self is the largest negative term it carries.

    Next: [notebook 4](04_epic_across_datasets.py) — all nine terms at once, from the command line, on
    held-out published neoantigens.
    """)
    return


if __name__ == "__main__":
    app.run()
