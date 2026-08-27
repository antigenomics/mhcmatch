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
    # 8 — From a scored table to a cassette (`mhcmatch.rank`, `mhcmatch.vector`)

    **What this demonstrates.** The two ends of the applied pipeline, on a **mock table** so it runs
    anywhere in seconds and every number is traceable to an input you can see. Notebooks 1–7 build
    the evidence; this one spends it.

    ```
    scored table  ->  rank  ->  units  ->  screen  ->  select  ->  order  ->  epitope_map
                     (what      (long     (safety)    (how many   (spacer +   (the map a
                      to keep)   windows)              per        ordering)    viewer draws)
                                                       allotype)
    ```

    **What you should conclude.**

    1. **`rank` and `vector` do not consume the same thing.** `rank` emits *minimal epitopes*; a
       cassette unit is the **long window around the mutation**. Feeding a ranking straight into
       `vector` is the single most common way to build a cassette with no flanks to process.
    2. **`screen` excludes, it never down-ranks.** Capacity spent on a unit that has to be
       withdrawn is capacity not spent on a safe one, so it runs *before* `select`.
    3. **The spacer is chosen, not assumed** — `order` tries **no spacer first** and picks what
       minimises predicted junctional binding.
    4. **Junctions manufacture epitopes that are in no unit.** `epitope_map` marks them `unit=0`;
       they are an artefact of assembly and belong to no gene.

    > Class I / CD8 throughout. The safety screen and the near-exact known-antigen columns are
    > scoped to class I — see :doc:`safety`.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.1 A mock scored table

    `rank_table` reads the pipeline `.scored.csv` schema — `epitope`, `best_allele`, `tpm`,
    `gene_name`, and whatever `score` the upstream tool produced. Six rows, written to a temp file,
    so nothing here depends on a cohort you do not have.
    """)
    return


@app.cell
def _():
    import csv
    import os
    import tempfile

    MOCK_ROWS = [
        # epitope,       best_allele,    tpm,    gene_name, upstream score
        ("KLIFLDSRV",   "HLA-A*02:01",  "142.0", "SIM2",   "0.81"),
        ("APRTLVYLL",   "HLA-B*07:02",  "61.3",  "ERG",    "0.49"),
        ("FYSLQEVMT",   "HLA-A*24:02",  "31.7",  "TP53",   "0.62"),
        ("YMDGTMSQV",   "HLA-A*02:01",  "18.9",  "TYR",    "0.58"),
        ("RPHERNGFTVL", "HLA-B*07:02",  "0.4",   "MAGEA1", "0.31"),
        ("AAAWYLWEV",   "HLA-A*02:01",  "77.2",  "MOCK1",  "0.44"),
    ]

    _d = tempfile.mkdtemp(prefix="mhcmatch-nb08-")
    TABLE = os.path.join(_d, "mock.scored.csv")
    with open(TABLE, "w", newline="") as _fh:
        _w = csv.writer(_fh)
        _w.writerow(["epitope", "best_allele", "tpm", "gene_name", "score"])
        _w.writerows(MOCK_ROWS)

    print(f"{TABLE}\n")
    print(open(TABLE).read())
    return MOCK_ROWS, TABLE, os


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.2 `rank` — recompute rather than re-sort

    Given a `store`, presentation is **recomputed** with this package's own binder score, so the
    ordering is ours and not a re-sort of someone else's column. The upstream number is kept in
    `components["score_builtin"]` precisely so the two can be compared rather than conflated.

    `tumor=` is worth passing, and it now does two things. It picks the matched normal tissue
    `expr_norm` reads, and it sets the abundance floor `c` that **both** expression terms divide by.
    A tumour's floor sits at roughly half its matched normal's — SKCM 0.1600 against skin 0.3050 TPM
    — so leaving it out is not a neutral choice: the pooled floor puts `expr_lvl` about a unit low.
    Where the origin arrives as free text, `expression.resolve_context` maps it, and an
    unrecognised string raises rather than quietly returning a number from the wrong distribution.
    """)
    return


@app.cell
def _(TABLE):
    import mhcmatch
    from mhcmatch import rank as R

    from mhcmatch import mimicry as MM

    store = mhcmatch.Store.from_pmhc(tier="shortlist", species="human", classes=("mhc1",))

    # The fitted aggregate scores on the features it declares, so the three corpus channels have
    # to be supplied. The CLI builds this callable in `cli._aggregate_channels`; the geometry
    # (k, face mask, substitution kernel) comes off the artifact, not a module default.
    def channels(peptides):
        g = MM.corpus_geometry()
        spec = MM.corpus_spectrum(cls="mhc1", components=("thymus", "self", "viral"),
                                  k=g["k"], self_species="human", mask=g["mask"],
                                  kernel=g["kernel"])
        rows = MM.corpus_R(list(peptides), spec, cls="mhc1")
        return {f"C_corpus_{c}": [r.get(c, float("nan")) for r in rows]
                for c in ("thymus", "self", "viral")}

    ranked = R.rank_table(TABLE, store=store, tumor="BRCA", channels=channels)

    print(f"{'peptide':<12}{'allele':<14}{'gene':<8}{'score':>8}{'pres':>8}{'expr':>9}  known")
    for r in ranked:
        print(f"{r.peptide:<12}{r.allele:<14}{r.gene:<8}{r.score:>8.4f}{r.presentation:>8.3f}"
              f"{r.expression:>9.1f}  {r.known_epitope or '-'}")
    return R, mhcmatch, ranked, store


@app.cell
def _(mo):
    mo.md(r"""
    Two things to read off that table rather than the score column.

    - **`known_epitope` is an annotation, never a fitted term.** An exact match to a tested
      neoantigen or a known viral epitope is direct evidence, and burying it inside a weighted sum
      lets a mediocre model score dilute the one hard fact in the row. `rank` flags it beside the
      score.
    - **The ordering is not the upstream ordering.** Compare `score` against
      `components["score_builtin"]` below — the inputs were handed to us pre-sorted by a different
      tool.
    """)
    return


@app.cell
def _(ranked):
    print(f"{'peptide':<12}{'ours':>9}{'upstream':>10}")
    for _r in ranked:
        _b = _r.components.get("score_builtin")
        print(f"{_r.peptide:<12}{_r.score:>9.4f}{('-' if _b is None else f'{_b:.4f}'):>10}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.3 Units — the long window, not the minimal epitope

    **This is the join that goes wrong most often.** `rank` emitted 9-mers; a cassette unit is the
    ~27-mer around the mutation, because a long peptide has to be *processed* by the recipient's own
    machinery and a bare 9-mer has no flanks to be cut out of. Kissick et al.
    ([PMID 24690990](https://pubmed.ncbi.nlm.nih.gov/24690990/)) is the case in point: one 27-mer
    whose overlapping class-I/class-II pair replaced an exogenous helper epitope outright.

    In a real run `vector.units_from_context(rows, records)` builds these by joining the ranking to
    the context FASTA — the mutation position lives in the FASTA header, not in the ranking, so the
    two are combined rather than either being guessed at. Here we write the windows out literally.
    """)
    return


@app.cell
def _():
    from mhcmatch import vector as V

    UNITS = [
        #        27-mer window                          mut  gene      allotype        p
        V.Unit("SVSTSGDLKLIFLDSRVTEVTGYSFRP", 13, "SIM2",   "HLA-A*02:01", 0.81),
        V.Unit("GHTRAPRTLVYLLDKDGNSVFVQAGET", 13, "ERG",    "HLA-B*07:02", 0.49),
        V.Unit("QLLDMKAFYSLQEVMTNQNRWKGVPLQ", 13, "TP53",   "HLA-A*24:02", 0.62),
        V.Unit("EYVIKVSARYMDGTMSQVQGSAKQRLL", 13, "TYR",    "HLA-A*02:01", 0.58),
        V.Unit("ARNDCQEGHILKMFPSTWYVAKLMNPQ", 13, "MOCK1",  "HLA-A*02:01", 0.44),
    ]
    for _u in UNITS:
        print(f"{_u.gene:<8}{_u.allele:<14}p={_u.p:.2f}  {_u.peptide}")
    return UNITS, V


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.4 `screen` — exclusion, before capacity is spent

    `risk(unit, registers) -> [reason, ...]`, empty meaning safe. It is **injected**, which is why
    the layout logic needs no proteome index; in production you pass
    `vector.self_origin_risk(...)`, built from the human proteome and reference expression, and it
    asks two questions: is the unit's *own* target gene transcribed in a tissue that must not be
    attacked, and does any register coincide with a self peptide from an **unrelated** essential
    gene.

    Both questions come from events, not hypotheticals: an affinity-enhanced MAGE-A3 TCR killed two
    patients through a titin off-target invisible to binding prediction
    ([PMID 23770775](https://pubmed.ncbi.nlm.nih.gov/23770775/)), and a MAGE-A3/A9/A12 TCR caused
    two more deaths because MAGE-A12 is transcribed in brain
    ([PMID 23377668](https://pubmed.ncbi.nlm.nih.gov/23377668/)).

    Here a mock `risk` stands in, firing on one unit so the withdrawal is visible.
    """)
    return


@app.cell
def _(UNITS, V):
    def mock_risk(unit, registers):
        """Stand-in for vector.self_origin_risk: pretend TP53's window hits a heart-expressed gene."""
        if unit.gene == "TP53":
            return [{"reason": "unrelated self origin", "gene": "MOCK-TITIN",
                     "tissue": "Heart - Left Ventricle", "register": registers[0]}]
        return []

    kept, rejected = V.screen(UNITS, mock_risk)
    print(f"kept     {[u.gene for u in kept]}")
    for _u, _reg, _why in rejected:
        print(f"withdrew {_u.gene}: {_why['reason']} -- {_why['gene']} in {_why['tissue']} "
              f"(register {_reg})")
    return kept, mock_risk, rejected


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.5 `select` — capacity, and why `n0` has no default

    `select` grows each allotype while the next candidate still beats that allotype's own expected
    yield per slot. Diversification across allotypes then **follows from the arithmetic** rather than
    from a quota someone imposed.

    `n0` is per-allotype capacity and is **required with no default**, in the library and in the
    Nextflow module alike: nothing in the public record fits it, so the value is yours to defend —
    and it is recorded in the output so a reader knows what was assumed.
    """)
    return


@app.cell
def _(V, kept):
    for _n0 in (1.0, 2.0, 4.0):
        _sel = V.select(kept, n0=_n0)
        print(f"n0={_n0:<5} -> {len(_sel.units)} units  {[(u.gene, u.allele.split('-')[1]) for u in _sel.units]}")

    sel = V.select(kept, n0=2.0)
    return (sel,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.6 `order` — the spacer is a result, not a convention

    `order` scores every register spanning every junction against the recipient's own allotypes and
    picks the spacer and ordering minimising predicted junctional binding — trying **no spacer
    first**, because a spacer that buys nothing is sequence that still has to be manufactured and
    translated.

    `objective="sum"` totals the strongest binder at each junction (pVACvector's logic,
    [PMID 31907209](https://pubmed.ncbi.nlm.nih.gov/31907209/)). It has a real bias toward the
    shortest spacer, which is a property of the objective and not a bug — `"rate"` is
    length-neutral and needs a `binder_threshold`.
    """)
    return


@app.cell
def _(V, sel, store):
    binder = V.store_binder(store, ["HLA-A*02:01", "HLA-B*07:02"], cls="mhc1")
    cassette = V.order(sel.units, binder)

    print(f"spacer   {cassette.spacer!r}")
    print(f"order    {[u.gene for u in cassette.units]}")
    print(f"length   {len(cassette.sequence)} aa   junction cost {cassette.cost:.4f}")
    print(f"\n{cassette.sequence}")
    return binder, cassette


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.7 `epitope_map` — the artefact a viewer can draw

    One row per **unit**, **linker** and predicted **epitope**, in 1-based inclusive coordinates over
    `cassette.sequence`. Three properties worth naming:

    - **A heterozygote is duplicated by construction** — one row per *(peptide, allele)*, so an
      epitope presented by two of the recipient's allotypes is two rows, which is what a coverage
      count needs.
    - **Junction-spanning epitopes carry `unit=0`** and no gene. They exist only because of the
      assembly, and telling them apart from real ones is the whole point of mapping.
    - **Class II gets `core_start`/`core_end`** — the register-anchored 9-mer inside the ligand —
      and `overlaps` links a class-II row to the class-I rows inside it, which is how `self_help`
      per unit is decided.
    """)
    return


@app.cell
def _(V, cassette, store):
    ranker1 = V.store_ranker(store, ["HLA-A*02:01", "HLA-B*07:02"], cls="mhc1")
    features = V.epitope_map(cassette, ranker1=ranker1, threshold=2.0)

    print(f"{'id':<5}{'kind':<9}{'start':>6}{'end':>5}  {'cls':<6}{'allele':<14}{'gene':<8}{'unit':>5}  seq")
    for f in features:
        print(f"{f.id:<5}{f.kind:<9}{f.start:>6}{f.end:>5}  {f.cls or '-':<6}{f.allele or '-':<14}"
              f"{f.gene or '-':<8}{f.unit:>5}  {f.seq}")
    return features, ranker1


@app.cell
def _(features):
    _junctional = [f for f in features if f.kind == "epitope" and f.unit == 0]
    print(f"{len(_junctional)} junction-spanning epitope(s) -- in no unit, belonging to no gene:")
    for _f in _junctional:
        print(f"  {_f.id}  {_f.seq}  {_f.allele}  %rank {_f.rank:.2f}  at {_f.start}-{_f.end}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.8 Writing it out

    `write_map` emits the flat TSV and a JSON a viewer can draw from **without recomputing
    anything** — the summary counts travel with the rows.
    """)
    return


@app.cell
def _(V, cassette, features, os):
    import json
    import tempfile as _tf

    _out = _tf.mkdtemp(prefix="mhcmatch-nb08-map-")
    _tsv = os.path.join(_out, "cassette.map.tsv")
    _json = os.path.join(_out, "cassette.map.json")
    summary = V.write_map(cassette, features, tsv_path=_tsv, json_path=_json)

    print(json.dumps(summary, indent=2)[:900])
    print(f"\n-- {_tsv}")
    print("\n".join(open(_tsv).read().splitlines()[:6]))
    return json, summary


@app.cell
def _(mo):
    mo.md(r"""
    ## 8.9 Back to nucleotides

    `back_translate` is **not** a codon optimiser. It fixes the two failure modes specific to a
    *concatemer* — the m1-pseudouridine +1-frameshift slippery motif (`slippery_sites`, whose seams
    the designer chose) and synthesis-hostile homopolymers, which spacers like `AAA` manufacture
    directly — and leaves GC content, secondary structure and CpG to a manufacturer's own tooling.
    `translate` exists so "synonymous" stays checkable rather than asserted.
    """)
    return


@app.cell
def _(V, cassette):
    cds = V.back_translate(cassette.sequence)
    print(f"{len(cds)} nt for {len(cassette.sequence)} aa")
    print(f"round-trips: {V.translate(cds) == cassette.sequence}")
    print(f"slippery sites remaining: {V.slippery_sites(cds)}")
    print(f"\n{cds[:120]}...")
    return (cds,)


@app.cell
def _(mo):
    mo.md(r"""
    ## What to carry away

    | step | the decision it makes | the trap it avoids |
    |---|---|---|
    | `rank_table` | recompute presentation, keep the upstream score beside it | re-sorting someone else's column and calling it a ranking |
    | units | the ~27-mer window, not the 9-mer | a cassette of minimal epitopes with no flanks to process |
    | `screen` | withdraw on safety **first** | spending capacity on a unit that gets withdrawn |
    | `select` | per-allotype saturation, `n0` explicit | a diversity quota nobody can defend |
    | `order` | spacer and order by junctional binding, **no spacer tried first** | a conventional `GPGPG` that buys nothing |
    | `epitope_map` | one row per (peptide, allele), junctions marked | counting a heterozygote once, or a junction artefact as a hit |

    Real inputs replace exactly two things: `units_from_context(rows, records)` instead of the
    hand-written `Unit`s, and `vector.self_origin_risk(...)` instead of `mock_risk`. Everything
    else is unchanged — which is the point of injecting `risk` and `binder`.
    """)
    return


if __name__ == "__main__":
    app.run()
