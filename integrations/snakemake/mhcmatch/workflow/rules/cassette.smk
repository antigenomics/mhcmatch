"""Select the units, assemble the construct, and score the whole cohort at once."""


def _alleles_in(w):
    if config.get("alleles"):
        return []
    return [f"{OUT}/alleles/{w.sample}.mhc1.txt"] if "typing" in SAMPLES_MAP[w.sample] else []


def _universe(w, input):
    """`--universe` is the DENOMINATOR coverage is reported against.

    Without it, coverage is taken over the labels the cassette happens to carry and cannot see the
    allotype it missed entirely -- which is the number a designer is asking for.
    """
    if config.get("alleles"):
        return f"--universe '{config['alleles']}'"
    return f"--universe \"$(cat {input.alleles[0]})\"" if input.alleles else ""


def _unit_source(w, input):
    """Where the assembly step reads each unit's LONG window from: `--context`, or a column.

    **There is no third option, and the fallback is not one.** A cassette unit is the ~27 aa window
    around the variant. `mhcmatch cassette`'s own `--unit-column` fallback is `peptide`, which on a
    reranked table is the MINIMAL epitope -- a 9-mer loads onto any cell without costimulation and
    is the TOLERISING configuration. So with neither a `windows` FASTA nor a configured
    `vector.unit_column` this stops the run and names the sample, rather than quietly building a
    construct out of minimal epitopes.
    """
    if input.context:
        return f"--context {input.context[0]}"
    col = config.get("vector", {}).get("unit_column")
    if not col:
        raise WorkflowError(
            f"sample {w.sample!r} ({w.arm} arm): no `windows` FASTA to pass as `--context` and no "
            "`vector.unit_column` set, so nothing names the long window a cassette unit has to be. "
            "The fallback would be `peptide` -- the MINIMAL epitope, which is the tolerising "
            "configuration -- so this is refused rather than defaulted. Give the sample a `windows` "
            "path in the samplesheet, or set `vector.unit_column` to the column carrying the ~27 aa "
            "window (`context_peptide` in integrations/fixtures/).")
    return f"--unit-column {col}"


rule mhcmatch_cassette_select:
    """Choose k epitopes by maximising the certainty-equivalent objective.

    Pass the WHOLE pool, not a shortlist: binding and expression carry the two largest coefficients
    in the model, so a pool already cut on them has no range left along the axes being traded.
    """
    input:
        pool=lambda w: pool_for(w.arm, w.sample),
        alleles=_alleles_in,
    output:
        f"{OUT}/{{arm}}/{{sample}}.vaccine.units.tsv",
    params:
        k=lambda w: config["cassette"]["k"],
        tol=lambda w: opt("--tol", config["cassette"]["tol"]),
        universe=_universe,
        scol=lambda w: opt("--score-column", score_column(w.arm)),
        block=lambda w: opt("--block-live", config["cassette"]["block_live"]),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch cassette select --candidates {input.pool} -k {params.k} {params.tol} "
        "{params.universe} {params.scol} {params.block} --passthrough --out {output}"


rule mhcmatch_cassette:
    """Order the chosen units, screen them, choose the linker, back-translate, emit the map."""
    input:
        units=f"{OUT}/{{arm}}/{{sample}}.vaccine.units.tsv",
        # The window FASTA on EITHER arm when the samplesheet gives one. A unit is the long window
        # whichever arm scored it, and `units_from_context` joins on the minimal epitope the units
        # table carries, so there is no reason for the two arms to read a different unit.
        context=lambda w: [SAMPLES_MAP[w.sample]["fasta_mhc1"]]
        if "fasta_mhc1" in SAMPLES_MAP[w.sample] else [],
        alleles=_alleles_in,
    output:
        report=f"{OUT}/{{arm}}/{{sample}}.cassette.tsv",
        protein=f"{OUT}/{{arm}}/{{sample}}.cassette.faa",
        cds=f"{OUT}/{{arm}}/{{sample}}.cassette.fna",
    params:
        # `order` and not `build`: `cassette select` has already chosen exactly k units and `build`
        # would re-select them under `--n0`, which is a different question. `order` keeps the
        # junction sweep, the back-translation and the safety screen.
        verb="order",
        alleles=lambda w, input: (f"--alleles '{config['alleles']}'" if config.get("alleles")
                                  else (f"--alleles \"$(cat {input.alleles[0]})\""
                                        if input.alleles else "")),
        ctx=_unit_source,
        screen=lambda w: flag("--screen", config["vector"]["screen"]),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch cassette {params.verb} --candidates {input.units} {params.ctx} "
        "{params.alleles} {params.screen} --fasta {output.protein} --fasta-nt {output.cds} "
        "--out {output.report}"


# ---------------------------------------------------------------------------------------------
# The cohort step. This is the one rule whose shape matters more than its command.
# ---------------------------------------------------------------------------------------------
def _all_units(w):
    return expand(f"{OUT}/{{arm}}/{{sample}}.vaccine.units.tsv", arm=w.arm, sample=SAMPLES)


def _all_pools(w):
    return [pool_for(w.arm, s) for s in SAMPLES]


rule mhcmatch_cassette_score:
    """One offset over EVERY donor in the run, and that is the whole point of the rule.

    `rank` anchors `p_response` on the batch it is handed, so a per-donor fit makes every donor's
    mean candidate probability equal the declared prevalence whatever their pool holds -- measured
    on 7,261 TCGA donors, every pool mean lands on 0.060163 with a standard deviation of 2.75e-17.
    Two donors' numbers are then the same number and a cross-donor triage built on them reads noise.

    **`expand()` over the module-level SAMPLES, never a glob over the output directory.** A glob is
    evaluated against whatever exists when the DAG is built, so it would silently fit the offset over
    a SUBSET -- reintroducing exactly the defect above while looking like it had been fixed. `expand`
    makes this rule unable to start until every donor's units exist.

    No `checkpoint`: the sample set is known from the input directory before the DAG is built, so a
    checkpoint would add re-evaluation machinery for nothing.
    """
    input:
        units=_all_units,
        pools=_all_pools,
    output:
        f"{OUT}/{{arm}}/cohort.cassette_score.tsv",
    params:
        scol=lambda w: opt("--score-column", score_column(w.arm)),
        prev=lambda w: opt("--prevalence", config.get("prevalence")),
        rho=lambda w: opt("--rho", config["cassette"]["rho"]),
        per=lambda w: flag("--per-donor-offset", config["cassette"]["per_donor_offset"]),
        block=lambda w: opt("--block-live", config["cassette"]["block_live"]),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch cassette score --cassettes {input.units} --pool {input.pools} "
        "{params.scol} {params.prev} {params.rho} {params.per} {params.block} --out {output}"


rule mhcmatch_cassette_report:
    """One self-contained HTML page over the same cohort. No jinja2, no matplotlib, no CDN."""
    input:
        units=_all_units,
        pools=_all_pools,
    output:
        f"{OUT}/{{arm}}/cohort.cassette_report.html",
    params:
        scol=lambda w: opt("--score-column", score_column(w.arm)),
        prev=lambda w: opt("--prevalence", config.get("prevalence")),
    threads: 1
    resources:
        mem_mb=2000,
        runtime=20,
    conda:
        "../../envs/mhcmatch.yaml"
    shell:
        "mhcmatch cassette report --cassettes {input.units} --pool {input.pools} "
        "{params.scol} {params.prev} --out {output}"
