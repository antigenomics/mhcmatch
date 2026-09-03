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
    # 11 — Linkers, and the mRNA (`mhcmatch.vector`)

    **What this demonstrates.** The last step of the cassette pipeline: what to put *between* the
    units, and how the amino-acid design becomes a molecule. Notebook 8 selects and orders units;
    this one picks the linker and builds the mRNA around them. Nothing here downloads anything or
    predicts anything — the payload is four hand-written units, so every number is checkable.

    **What you should conclude.**

    1. A linker is a **named preset with provenance**, not a string. `LINKERS` carries the family,
       the class each is intended for, and where it comes from — including that the GS 10-mer is a
       reconstruction of a format that was described rather than published.
    2. The preset table deliberately **does not rank itself**. The two published mechanisms that
       would settle a class-I ranking act at different positions and point opposite ways, so which
       linker wins is a measurement against the recipient's own allotypes, not a citation.
    3. The linker is not free at the nucleotide level either. Back-translating the *whole* reading
       frame in one pass is what lets homopolymer avoidance act on the seams the linker created —
       and the seams are where a concatemer's problems are.
    4. `mrna()` returns a parts map that **tiles the molecule exactly**, so an element that went
       missing is a gap and one placed twice is an overlap. A length check catches neither.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 1. The preset table

    Six families, each buying something different. `cls` is the class a linker is *intended* for,
    which is provenance and not a measurement.
    """)
    return


@app.cell
def _():
    from mhcmatch import vector as V

    for _name, _L in V.LINKERS.items():
        print(f"{_name:<8} {_L.sequence or '-':<11} {len(_L):>2}aa  "
              f"{_L.family:<19} {_L.cls:<5} {_L.note[:58]}")
    return (V,)


@app.cell
def _(mo):
    mo.md(r"""
    ### Why it does not rank itself

    `GPGPG` is the only linker with causal evidence behind it, and the evidence is **class II**:
    unspaced concatenation of four HLA-DR epitopes created a junctional epitope that suppressed the
    response to all four, and inserting `GPGPG` restored all four (Livingston *et al.*,
    *J Immunol* 2002, PMID 12023344). Against that, all six orderings of three Fel d I regions
    produced no detectable junctional response at all (Rogers *et al.*, *Mol Immunol* 1994,
    PMID 7521933). Both results are real.

    For class I the two candidate mechanisms act **at different positions**:

    | argues | for | mechanism |
    |---|---|---|
    | Martin-Galiano & Lopez, *PLoS One* 2019 (PMID 30645615) | Gly/Pro | glycine and proline are abundant in the C-terminal regions class-I ligands are cleaved from, which is a processing argument for them |
    | Bergmann *et al.*, *J Immunol* 1996 (PMID 8871618) | against Gly/Pro | glycine and proline flanking a class-I epitope inhibit recognition of the epitope on their **amino-terminal** side; the ratio of two responses from one construct moved up to **50-fold** with flanking context |

    A preset chosen on either citation alone is a linker chosen by reputation. `order()` scores each
    candidate against the recipient's own allotypes instead, and *that* is what selects — pinning
    with `linker=` is for when the construct format is already decided and only the layout is open.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 2. Assembling, with the format already decided

    `assemble()` lays units out in the order given, joined by a named preset, and predicts nothing.
    `boundaries` tiles the units, so the gaps between them are exactly the linkers.
    """)
    return


@app.cell
def _(V):
    units = [V.Unit("MTEYKLVVVGAKVAELVHFLGDPTIE", 13, "KRAS", "HLA-A*02:01", 0.31),
             V.Unit("KLQEEIPVLSIINFEKLAKQVWRTAY", 13, "OVA", "HLA-A*02:01", 0.42),
             V.Unit("QPRVLTQEQAGTLSHFWDNPYTKQRA", 13, "TP53", "HLA-B*07:02", 0.27),
             V.Unit("EDLLKYYSQLNPRTGSWDMQNLKAAV", 13, "PIK3CA", "HLA-B*07:02", 0.19)]

    for _name in ("none", "AAY", "GS10", "GPGPG"):
        _c = V.assemble(units, _name)
        _gaps = {_c.sequence[a[1]:b[0]]
                 for a, b in zip(_c.boundaries, _c.boundaries[1:])}
        print(f"{_name:<7} {len(_c.sequence):>4} aa   linker={_c.spacer or '-':<11} "
              f"gaps={_gaps}")
    return (units,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 3. The molecule

    `mrna()` builds the construct and returns what it built. The backbone below is a **placeholder**
    — the library supplies no UTR, no signal peptide and no tail, because those belong to a
    particular vector and a plausible invented one is worse than none.
    """)
    return


@app.cell
def _(V, units):
    # Every backbone element below is a DUMMY, chosen so it is obviously not a real sequence. The
    # library supplies none of them, and an example that pasted in a plausible-looking signal
    # peptide would be read as a recommendation. Substitute your own vector's elements.
    m = V.mrna(units, "GS10",
               leader="M" + "A" * 18,                 # dummy signal peptide, 19 aa
               trailer="H" * 24,                      # dummy trafficking domain, 24 aa
               utr5="ACGT" * 8, utr3="ACGT" * 4, poly_a=100)

    for _k, _v in m.checks.items():
        print(f"{_k:>20}  {_v}")
    print()
    for _p in m.parts:
        print(f"{_p['kind']:<8} {_p['name'][:22]:<22} {_p['start'] + 1:>5}-{_p['end']:<5} "
              f"{_p['end'] - _p['start']:>4} nt")
    return (m,)


@app.cell
def _(mo):
    mo.md(r"""
    Two properties are worth reading off that table rather than trusting.

    **The parts tile the molecule.** Consecutive, non-overlapping, and the concatenation of every
    slice is the whole sequence. That is what makes the record auditable: a missing element shows up
    as a gap and a duplicated one as an overlap, and a length check sees neither.

    **`translates` is the check that must hold.** The coding sequence, read back in the frame the
    assembled construct sets it in, gives exactly `protein`. It is what catches a frame broken by a
    supplied 5' element — the failure mode that is otherwise invisible.
    """)
    return


@app.cell
def _(m):
    assert "".join(m.slice_of(p) for p in m.parts) == m.sequence
    assert all(a["end"] == b["start"] for a, b in zip(m.parts, m.parts[1:]))
    print(f"{len(m.parts)} parts tile {len(m.sequence)} nt exactly; "
          f"translates={m.checks['translates']}")
    print(m.as_rna()[:96], "...")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## 4. The linker is not free in nucleotides either

    `AAA` is three alanines, whose codons all begin `GC`; a lysine either side of the seam is
    `AAG`/`AAA`. A run like that is a property of the **junction**, not of either unit, so a
    per-unit back-translation cannot see it — which is why the whole reading frame is
    back-translated in one pass. The column below is what that buys: every linker comes back at the
    `MAX_HOMOPOLYMER` target of 4 and with no m1Ψ +1-frameshift motif left, *including* the ones
    whose own residues manufacture the run. What the linker does move is length and GC, and GC is
    reported rather than managed.
    """)
    return


@app.cell
def _(V, units):
    print(f"{'linker':<8} {'nt':>5} {'GC':>7} {'max run':>8} {'slippery':>9}")
    for _name in ("none", "G", "AAA", "AAY", "GS10", "G4S2", "GPGPG", "EAAAK"):
        _m = V.mrna(units, _name)
        _c = _m.checks
        print(f"{_name:<8} {_c['cds_nt']:>5} {_c['gc']:>7.3f} "
              f"{_c['longest_homopolymer']:>8} {_c['slippery_sites']:>9}")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## What this is not

    Not a codon optimiser. It fixes the two things that make a *polyepitope* construct fail where a
    natural ORF would not, and leaves GC content, secondary structure, splice sites and CpG to the
    manufacturer's tooling — the GC column above is reported, not managed.

    Not a recommendation of a linker. The panel is scored, and which one wins for a given recipient
    is `order()`'s answer on that recipient's allotypes.

    ### Where the numbers come from

    Nothing here. The four units are written by hand so every figure is reproducible without a
    download; the backbone elements are placeholders. Measured cassette figures live in
    [`2026-mhcmatch-code`](https://github.com/repseq/2026-mhcmatch-code) (private; released to reviewers) under
    `bench/results/`.
    """)
    return


if __name__ == "__main__":
    app.run()
