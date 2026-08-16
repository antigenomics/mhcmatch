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
    # 3 — T-cell precursor frequency (`mhcmatch.precursor`)

    **What this demonstrates.** How much naive-repertoire mass can see a given epitope. The estimand is

    $$F(e) = \sum_{\tau \in C_e} \pi(\tau)$$

    — the probability that a random naive-repertoire junction falls in the cognate set $C_e$. This is
    the continuous quantity behind the word "immunogenic".

    `mhcmatch.precursor` offers **three estimators of the same $F(e)$**, in increasing order of model
    commitment, and the notebook runs all three on one epitope's real TCRs:

    | estimator | what it sums over | bias |
    |---|---|---|
    | `observed_mass` | the junctions actually recorded | a **strict lower bound**, and a biased one |
    | `ball_mass` | the **union** of Hamming-$r$ balls around them | closer, but still anchored on the observed sample |
    | `motif_mass` | every junction matching a degenerate motif | scores a **set**, so no coverage bias at all |

    **What you should conclude.**

    1. `observed_mass` is a floor, not an estimate. VDJdb samples cognate TCRs **size-biased by Pgen**
       — a TCR enters the record roughly in proportion to its repertoire frequency — so the recorded
       members are systematically the high-Pgen ones and the deficit does not shrink with more studies
       at the same depth.
    2. `ball_mass` returns a **union**, not a sum, and reports the `overlap` you would have invented by
       summing. Section 3.4 shows that overlap climbing as the specificity group gets denser.
    3. **CDR3 is not junction.** An anchor-stripped IMGT CDR3 scores exactly `0.0` with no error.
       `check_junctions` exists so that failure is loud instead of silent.

    Requires the optional `[precursor]` extra (`pip install 'mhcmatch[precursor]'`), which pulls in
    `vdjtools` for the recombination model. Nothing here reimplements Pgen.

    Data: VDJdb, fetched from the public HuggingFace dataset `isalgo/airr_benchmark`.
    """)
    return


@app.cell
def _():
    import collections
    import csv
    import gzip
    import time

    from huggingface_hub import hf_hub_download

    from mhcmatch import precursor

    csv.field_size_limit(10**7)  # VDJdb slim has very wide method/meta fields

    VDJDB_REPO = "isalgo/airr_benchmark"
    VDJDB_FILE = "vdjdb/vdjdb-2026-06-11-ZENODO/vdjdb.slim.txt.gz"
    # Two epitopes whose VDJdb records are artefactual; excluded from any census.
    EXCLUDE = {"SLLMWITQV", "KLGGALQAK"}

    _t0 = time.time()
    vdjdb_path = hf_hub_download(repo_id=VDJDB_REPO, repo_type="dataset", filename=VDJDB_FILE)
    print(f"VDJdb fetched/cached in {time.time() - _t0:.1f} s")

    _t0 = time.time()
    trb_rows = []
    with gzip.open(vdjdb_path, "rt") as _fh:
        for _r in csv.DictReader(_fh, delimiter="\t"):
            if _r["species"] != "HomoSapiens" or _r["gene"] != "TRB":
                continue
            if _r["antigen.epitope"] in EXCLUDE:
                continue
            # NB: the column is named `cdr3` but holds JUNCTIONS -- see section 3.2.
            trb_rows.append(
                (_r["antigen.epitope"], _r["cdr3"], _r["v.segm"], _r["j.segm"], _r["mhc.a"])
            )
    print(f"parsed {len(trb_rows)} human TRB records in {time.time() - _t0:.1f} s")

    census = collections.Counter(_row[0] for _row in trb_rows)
    for _ep, _n in census.most_common(6):
        print(f"  {_ep:<12s} {_n:>6d} TRB records")
    return collections, precursor, time, trb_rows


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.1 One epitope's cognate set

    `GILGFVFTL` (influenza A M1 58-66, HLA-A\*02:01) — the most-studied public specificity there is.
    """)
    return


@app.cell
def _(trb_rows):
    EPITOPE = "GILGFVFTL"

    cognate = sorted({(_c, _v, _j) for _e, _c, _v, _j, _m in trb_rows if _e == EPITOPE})
    all_junctions = sorted({_c for _c, _v, _j in cognate})
    print(f"{EPITOPE}: {len(cognate)} unique (junction, V, J) records")
    print(f"           {len(all_junctions)} unique junctions")
    print(f"example:   {cognate[0]}")
    return EPITOPE, all_junctions, cognate


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.2 The CDR3-vs-junction trap

    **VDJdb's column is named `cdr3` but holds junctions** — Cys104 and Phe118/Trp118 included. An
    IMGT CDR3 is the same string with those two anchors stripped, so it is two residues shorter.

    Feed a stripped CDR3 to the recombination model and it does not raise: it returns **exactly
    `0.0`**, because no rearrangement produces a junction without its anchors. A mis-typed input
    therefore reports a precursor frequency of zero instead of failing. `check_junctions` splits the
    input on the conserved anchors so the caller can report the dropped count.
    """)
    return


@app.cell
def _(all_junctions, mo, precursor):
    model = precursor.load_model("TRB")  # bundled OLGA model; source="olga" is the default for null work

    junctions, suspect = precursor.check_junctions(all_junctions)
    print(f"check_junctions: {len(junctions)} junctions, {len(suspect)} suspect")
    print(f"suspect entries: {suspect}")

    _j = junctions[0]
    _stripped = _j[1:-1]  # what an IMGT CDR3 column would contain
    _mass_j = precursor.observed_mass(model, [_j])
    _mass_s = precursor.observed_mass(model, [_stripped])
    print()
    print(f"junction   {_j:<20s} Pgen = {_mass_j:.6e}")
    print(f"CDR3       {_stripped:<20s} Pgen = {_mass_s:.6e}   <- silent zero, no exception")

    mo.md(
        f"""
    | input | string | length | `observed_mass` |
    |---|---|--:|--:|
    | junction (what VDJdb's `cdr3` column holds) | `{_j}` | {len(_j)} | **{_mass_j:.3e}** |
    | IMGT CDR3 (anchors stripped) | `{_stripped}` | {len(_stripped)} | {_mass_s:.1f} |

    Always `check_junctions` first, and report the dropped count rather than letting it vanish. The
    {len(suspect)} suspect entries above are real VDJdb rows that do not end on the conserved
    Phe/Trp.
    """
    )
    return junctions, model


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.3 `observed_mass` — the lower bound

    A plain sum of Pgen over the recorded junctions. It is *strictly* a lower bound: any cognate TCR
    that was never sequenced contributes nothing.

    The bias is not random. VDJdb samples cognate TCRs size-biased by Pgen, so what got recorded is
    the high-Pgen tail of the cognate set, and the deficit does not shrink by adding more studies at
    the same sequencing depth.
    """)
    return


@app.cell
def _(EPITOPE, junctions, model, precursor, time):
    _t0 = time.time()
    observed = precursor.observed_mass(model, junctions)
    print(f"observed_mass over {len(junctions)} junctions: {observed:.6e}  ({time.time() - _t0:.1f} s)")
    print(f"  i.e. ~{observed * 1e6:.1f} in 10^6 naive TRB rearrangements are a *recorded* {EPITOPE} TCR")

    # The mass is concentrated: a handful of highly generatable junctions dominate the sum.
    _per_seq = sorted(((precursor.observed_mass(model, [_j]), _j) for _j in junctions), reverse=True)
    print(f"  top junction: {_per_seq[0][1]} at {_per_seq[0][0]:.3e}")
    print(f"  top 10 of {len(junctions)} junctions carry "
          f"{sum(_m for _m, _ in _per_seq[:10]) / observed:.1%} of the total mass")
    return (observed,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.4 `ball_mass` — the **union** of Hamming-1 balls, not the sum

    Cognate TCRs are near-duplicates by construction: that is what a specificity group *is*. So their
    Hamming-1 balls intersect, and adding per-sequence ball masses counts the intersections more than
    once. `ball_mass` enumerates the union (via `seqtree.neighbourhood_union`) and scores it once.

    It returns `naive_sum` alongside, and `overlap = 1 - union/naive_sum` — the fraction of mass that
    double-counting would have invented. That number is worth reporting rather than hiding: it is a
    direct measurement of how tight the specificity group is.

    Below, three nested subsets of the length-13 junctions, at increasing Hamming radius from one
    seed. As the group gets denser, the overlap climbs.

    V and J are deliberately **not** accepted by `ball_mass`: a substituted neighbour need not keep the
    centre's V/J assignment, so conditioning the ball on the centre's call would be wrong. It
    marginalises.
    """)
    return


@app.cell
def _(junctions, mo, model, precursor, time):
    SEED = "CASSIRSSYEQYF"


    def _hamming(a, b):
        return sum(x != y for x, y in zip(a, b)) if len(a) == len(b) else 1 << 30


    _l13 = [_s for _s in junctions if len(_s) == 13]
    print(f"length-13 junctions: {len(_l13)}")

    ball_rows = []
    for _r in (1, 2, 3):
        _grp = [_s for _s in _l13 if _hamming(SEED, _s) <= _r]
        _t0 = time.time()
        _b = precursor.ball_mass(model, _grp, r=1)
        print(
            f"  within Hamming {_r} of the seed: n={len(_grp):>3d}  "
            f"union={_b['union']:.3e}  overlap={_b['overlap']:.1%}  ({time.time() - _t0:.1f} s)"
        )
        ball_rows.append((_r, len(_grp), _b))

    _rows = "\n".join(
        f"| {_r} | {_n} | {_b['n_union']} | {_b['union']:.3e} | {_b['naive_sum']:.3e} | {_b['overlap']:.2%} |"
        for _r, _n, _b in ball_rows
    )
    mo.md(
        f"""
    **Seed `{SEED}`, length-13 cognate junctions**

    | Hamming radius from seed | n centres | n sequences in the union | `union` | `naive_sum` | `overlap` |
    |--:|--:|--:|--:|--:|--:|
    {_rows}

    At {ball_rows[-1][1]} centres the naive sum invents **{ball_rows[-1][2]["overlap"]:.1%}** of the
    mass. The union is the correct object; `naive_sum` is only there so you can see what you avoided.
    """
    )
    return (ball_rows,)


@app.cell
def _(mo):
    mo.md(r"""
    ## 3.5 `motif_mass` — degenerate Pgen, no coverage bias

    A VDJdb cluster PWM is V/J/length-pinned, so thresholding it per position gives exactly a
    **per-position set of permitted residues** — which is what `motif_mass` takes. One masked DP pass
    returns the summed Pgen of every junction the motif matches: no enumeration, no
    inclusion-exclusion, and no dependence on which cognate TCRs happened to be sequenced.

    First the semantics, on a single sequence.
    """)
    return


@app.cell
def _(junctions, model, precursor):
    _j = junctions[0]

    _pinned = precursor.motif_mass(model, list(_j))  # every position pinned == the single sequence
    _single = precursor.observed_mass(model, [_j])
    print(f"motif_mass(list({_j!r}))  = {_pinned:.6e}")
    print(f"observed_mass([{_j!r}])   = {_single:.6e}")
    print(f"identical: {abs(_pinned - _single) < 1e-30}")

    _wide = list(_j)
    _wide[5] = ""  # "" or "X" == wildcard at that position
    print()
    print(f"widening position 6 to a wildcard: {precursor.motif_mass(model, _wide):.6e}")
    print("(widening a position can only add mass)")
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### A real motif

    Take the V/J/length-pinned subgroup of the cognate set, and keep at each position every residue
    seen in at least 5% of its members. That is a thresholded cluster PWM, and `motif_mass` scores the
    whole set it describes in one call.
    """)
    return


@app.cell
def _(cognate, collections, mo, model, observed, precursor):
    MOTIF_V, MOTIF_J, MOTIF_L = "TRBV19*01", "TRBJ2-3*01", 13
    MOTIF_FLOOR = 0.05

    group = sorted({_c for _c, _v, _j in cognate if _v == MOTIF_V and _j == MOTIF_J and len(_c) == MOTIF_L})
    _cols = [collections.Counter(_s[_i] for _s in group) for _i in range(MOTIF_L)]
    allowed = []
    for _col in _cols:
        _tot = sum(_col.values())
        allowed.append("".join(sorted(_a for _a, _n in _col.items() if _n / _tot >= MOTIF_FLOOR)) or "X")

    group_observed = precursor.observed_mass(model, group)
    motif_marginal = precursor.motif_mass(model, allowed)
    motif_conditioned = precursor.motif_mass(model, allowed, v=MOTIF_V, j=MOTIF_J)
    group_ball = precursor.ball_mass(model, group, r=1)

    print(f"group: {MOTIF_V} / {MOTIF_J} / length {MOTIF_L}  ->  {len(group)} junctions")
    print(f"motif: {'-'.join(allowed)}")
    print(f"observed_mass(group)             = {group_observed:.6e}")
    print(f"ball_mass(group)['union']        = {group_ball['union']:.6e}")
    print(f"motif_mass(allowed)              = {motif_marginal:.6e}   (V/J marginalised)")
    print(f"motif_mass(allowed, v=..., j=...) = {motif_conditioned:.6e}   (V/J conditioned)")

    mo.md(
        f"""
    **Motif** (`{MOTIF_V}` / `{MOTIF_J}` / length {MOTIF_L}, {len(group)} observed members, residues at
    >= {MOTIF_FLOOR:.0%} of the column):

    `{"-".join(allowed)}`

    | estimator | value | scores over |
    |---|--:|---|
    | `observed_mass(group)` | {group_observed:.3e} | the {len(group)} recorded junctions |
    | `ball_mass(group)["union"]` | {group_ball["union"]:.3e} | their Hamming-1 union ({group_ball["n_union"]} sequences) |
    | `motif_mass(allowed)` | {motif_marginal:.3e} | every junction the motif matches, V/J marginalised |
    | `motif_mass(allowed, v=, j=)` | {motif_conditioned:.3e} | the same set, V/J conditioned |

    `motif_mass` and `ball_mass` estimate the same quantity by independent routes, so **their
    disagreement is itself an estimate of the missing mass**.

    **Marginal and conditioned masses are different quantities — never mix them in one comparison.**
    The marginal is the larger of the two; conditioning on the motif's own V/J is the right choice for
    a cluster motif (they are V/J-pinned by construction), while `ball_mass` has no choice but to
    marginalise.

    For scale: the whole recorded cognate set gives `observed_mass` = {observed:.3e}, and this single
    motif accounts for {motif_marginal:.3e} on its own.
    """
    )
    return


@app.cell
def _(EPITOPE, ball_rows, mo, observed):
    mo.md(
        f"""
    ---

    ## Summary for `{EPITOPE}`

    | estimator | value | interpretation |
    |---|--:|---|
    | `observed_mass` (all recorded junctions) | {observed:.3e} | strict lower bound, Pgen-size-biased |
    | `ball_mass` union (densest subset shown) | {ball_rows[-1][2]["union"]:.3e} | union of Hamming-1 balls; naive summing would have added {ball_rows[-1][2]["overlap"]:.1%} |

    ### Practical notes

    * `load_model(locus, source=, organism=)` — **do not use `source="learned"`**. Those models were
      EM-fit on ~2k clonotypes without a gene-usage pseudocount, so most bundled TRB V alleles have
      `P(V) = 0` and any junction using them scores 0. Mouse is available only under `source="arda"`,
      and only for TRA/TRB.
    * `v`/`j` arguments take **allele-resolution** names (`TRBV27*01`, not `TRBV27`).
    * `threads=0` lets the Pgen backend pick; the DP, the closed Hamming-1 ball and the degenerate DP
      all live in `vdjtools`, and the neighbourhood enumeration in `seqtree`.
    * `python -m mhcmatch.precursor` runs the module's own self-check (and skips cleanly without
      `vdjtools`).
    """
    )
    return


if __name__ == "__main__":
    app.run()
