// mhcmatch as nf-core-style local processes. Nine, in pipeline order:
//
//   MHCMATCH_ALLELES          an HLA typing file -> the allele list every other process takes
//   MHCMATCH_PREDICT          variant windows    -> per-allele presentation + affinity
//   MHCMATCH_RANK             windows            -> the fitted EPIC aggregate, one ordered table
//   MHCMATCH_RERANK           a caller's OWN table -> the same aggregate, appended to it
//   MHCMATCH_NEOAG            peptides           -> proximity to the tested-neoantigen database
//   MHCMATCH_MIMICRY          peptides           -> the signed self/viral/thymus channels
//   MHCMATCH_CASSETTE_SELECT  a pool             -> the k units to manufacture
//   MHCMATCH_CASSETTE         ranked units       -> a screened polyepitope cassette, aa and CDS
//   MHCMATCH_CASSETTE_SCORE   EVERY donor        -> one shared calibration
//
// `../overlay/` attaches these to a host pipeline; `pipeline.nf` runs them from a samplesheet.
// `container`, `conda` and `publishDir` live in nextflow.config, in one selector, not per process.
//
// Two conventions worth knowing before editing:
//
//  * **No stub types a header.** Every stub asks the installed library for its own schema, because
//    a hand-copied header drifts -- this module shipped an 18-column `scored.csv` stub against a
//    57-column real table once. `predict.SCORED_COLUMNS`, `predict.NATIVE_COLUMNS`,
//    `rank.columns()` and `mimicry.NEOAG_COLUMNS` are the sources of truth.
//  * **Species follows `params.genome`**, mapped in nextflow.config via ext.args, so there is no
//    extra parameter to configure.

//: Is a boolean parameter on?
//:
//: **`--flag false` on the command line arrives as the STRING "false", which is truthy in Groovy**,
//: so the plain `params.x ? '--flag' : ''` idiom passes the flag a user just tried to disable. It
//: cannot be fixed in nextflow.config either: a config statement is evaluated BEFORE Nextflow
//: applies `--param`, so a coercion written there is overwritten by the value it exists to coerce.
//: The direction that matters is the reverse one -- somebody who believes they enabled
//: `--mhcmatch_cassette_screen` and did not gets a cassette with no safety check and no error.
def isOn(v) {
    v != null && !(v.toString().toLowerCase() in ['false', '0', 'no', 'null', ''])
}

//: The whitelist flags. Two lists, because they make different claims about a row: a gene hit says
//: the gene is of interest, an epitope hit is evidence about the peptide itself. `--keep-epitopes
//: builtin` loads a pre-built seqtree index off disk (~1 ms) -- nothing here builds an index or
//: writes a cache, so N concurrent samples cannot race.
def keepArgs() {
    def out = ''
    if (params.mhcmatch_keep_genes)    { out += "--keep-genes '${params.mhcmatch_keep_genes}' " }
    if (params.mhcmatch_keep_epitopes) {
        out += "--keep-epitopes '${params.mhcmatch_keep_epitopes}' "
        // Passed THROUGH, not folded: an out-of-range radius is refused by name by the CLI rather
        // than silently narrowed here.
        if (isOn(params.mhcmatch_keep_mismatch)) { out += "--keep-mismatch ${params.mhcmatch_keep_mismatch} " }
    }
    out
}

//: The three `rank`/`rerank` flags that are plain on/off.
def rankExtras() {
    (isOn(params.mhcmatch_rank_extended) ? '--extended ' : '') +
    (isOn(params.mhcmatch_rank_annotate) ? '--annotate ' : '') +
    (isOn(params.mhcmatch_rank_core)     ? '--core '     : '')
}

//: `human` / `mouse`, for a stub that has to name the species the real run would.
def speciesOf() {
    params.genome == 'GRCm39' ? 'mouse' : 'human'
}

//: The output prefix for a process BELOW the arm fork.
//:
//: The two arms run the same five processes over the same donor, so `${meta.id}` alone is two
//: tasks publishing one filename. The arm rides in `meta` rather than in six `withName:` selectors
//: -- a selector spelled as a bare process name matches one arm and silently misses the other,
//: which for MHCMATCH_CASSETTE once meant 8 GB instead of 48 and an OOM kill hours in.
def armPrefix(meta, ext) {
    ext ?: (meta.arm ? "${meta.id}.${meta.arm}" : "${meta.id}")
}

//: Which column of a pool holds the aggregate.
//:
//: **Do not leave this to the CLI fallback on the rerank arm.** `_cassette_rows` falls back to
//: `score`, and a caller's candidate table HAS one -- theirs -- so the arm would select on their
//: ranking while looking as though it selected on ours.
def scoreColumn(meta) {
    params.mhcmatch_cassette_score_column ?:
        (meta?.arm == 'rerank' ? "${params.mhcmatch_rerank_prefix}score" : '')
}


process MHCMATCH_ALLELES {
    tag "${meta.id}:${cls}"
    label 'process_single'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    // **The step whose absence is silent.** Every HLA typer writes the G-group form
    // (`A*01:01:01G`), the pseudosequence tables are keyed at two fields, and `Store._allele_set`
    // drops what it cannot find WITHOUT A WORD -- so a run handed a raw typing file scores against
    // an empty panel and exits 0. The class-II half is worse: `DQA1*05:01` on its own is not a
    // molecule, so the two rows of a typing file have to be joined.
    input:
    tuple val(meta), path(typing), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.alleles.txt"), emit: alleles

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mhcmatch alleles ${typing} --cls ${cls} ${task.ext.args ?: ''} \\
        --out ${prefix}.${cls}.mhcmatch.alleles.txt
    """

    stub:
    """
    printf '' > ${task.ext.prefix ?: meta.id}.${cls}.mhcmatch.alleles.txt
    """
}


process MHCMATCH_PREDICT {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    input:
    tuple val(meta), path(fasta), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.scored.csv"), emit: scored
    // `native_tsv`, not `native`: `native` is a Groovy reserved word, and Nextflow 21.10.6 fails
    // the WHOLE module with "Unexpected input: '{'" pointing at the enclosing `process {`, which
    // reads as a corrupt file rather than one bad identifier.
    tuple val(meta), val(cls), path("*.mhcmatch.native.tsv"), emit: native_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mhcmatch predict ${fasta} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${params.mhcmatch_tier} \\
        --rank-threshold ${params.mhcmatch_rank_threshold} \\
        ${keepArgs()}${isOn(params.mhcmatch_predict_core) ? '--core ' : ''}${task.ext.args ?: ''} \\
        --scored-csv ${prefix}.${cls}.mhcmatch.scored.csv \\
        --native ${prefix}.${cls}.mhcmatch.native.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c "from mhcmatch.predict import SCORED_COLUMNS as C; print(','.join(C))" \\
        > ${prefix}.${cls}.mhcmatch.scored.csv
    python -c "from mhcmatch.predict import NATIVE_COLUMNS as C; print('\\t'.join(C))" \\
        > ${prefix}.${cls}.mhcmatch.native.tsv
    """
}


process MHCMATCH_RANK {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    // `rank` reads the known-epitope sets, the mimicry references and the expression tables on top
    // of the ligand panel. The image bakes them (`bootstrap --reference`); a bare `bootstrap` image
    // reaches for HuggingFace from the compute node instead.
    input:
    tuple val(meta), path(input), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.ranked.tsv"), emit: ranked

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The immunological MODE, not the input shape -- `mhcmatch_rank_mode` is the shape. `--tumor`
    // is dropped in pathogen mode because `rank` refuses it there (undefined without a host
    // transcript) and would exit non-zero on every task in the arm.
    def epi   = params.mhcmatch_rank_epitope
    def tumor = (epi == 'neoantigen' && params.mhcmatch_tumor) ? "--tumor ${params.mhcmatch_tumor} " : ''
    // `p_response`'s anchor is a PRIOR about this cohort's candidate list, not a model output, so
    // it is a pipeline parameter. Left unset, the CLI uses TESLA's 37 of 615.
    def prev  = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence} " : ''
    """
    mhcmatch rank ${params.mhcmatch_rank_mode} ${input} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${params.mhcmatch_tier} \\
        --rank-threshold ${params.mhcmatch_rank_threshold} \\
        --epitope ${epi} \\
        --score ${params.mhcmatch_rank_score} \\
        ${keepArgs()}${tumor}${prev}${rankExtras()}${task.ext.args ?: ''} \\
        --out ${prefix}.${cls}.mhcmatch.ranked.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // **`cls`, `species` and `mode` reach the stub too, or the stub is a different header.** A run
    // emits the columns it COMPUTED: pathogen mode drops the wild-type and expression blocks, and a
    // class-II fit declares no corpus channels.
    """
    python -c "from mhcmatch import rank; print('\\t'.join(rank.columns(\\
        extended=${isOn(params.mhcmatch_rank_extended) ? 'True' : 'False'}, \\
        annotate=${isOn(params.mhcmatch_rank_annotate) ? 'True' : 'False'}, \\
        core=${isOn(params.mhcmatch_rank_core) ? 'True' : 'False'}, \\
        score='${params.mhcmatch_rank_score}', cls='${cls}', species='${speciesOf()}', \\
        mode='${params.mhcmatch_rank_epitope}')))" > ${prefix}.${cls}.mhcmatch.ranked.tsv
    """
}


process MHCMATCH_RERANK {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    // `rank pairs --passthrough`: the caller's OWN table comes back with every column it arrived
    // with, in its own order, plus this model's under `--prefix`, re-ordered by the aggregate. Not
    // a join a caller can do afterwards -- `rank` splits a cell naming several alleles and the best
    // presenter stands for the row, so the output shares neither its length nor its allele column
    // with the input.
    //
    // `context` is the window FASTA the candidates were called on, and it is what makes
    // agretopicity and `d_occupancy` defined: a candidate table carries the mutant k-mer and
    // nothing the germline is recoverable from. Pass NO_FILE and every row is `wt_absent`.
    //
    // No `--alleles` here, and none is accepted: `rank pairs` refuses it, because the rows name
    // their own. No `--rank-threshold` either -- every row of the caller's table comes back, which
    // is the contract this arm exists for. The whitelists still apply; there they only ever set
    // `keep` / `keep_reason`.
    input:
    tuple val(meta), path(table), path(context), val(cls)

    output:
    tuple val(meta), val(cls), path("*.epitopes.mhcmatch.tsv"), emit: reranked

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def epi    = params.mhcmatch_rank_epitope
    def tumor  = (epi == 'neoantigen' && params.mhcmatch_tumor) ? "--tumor ${params.mhcmatch_tumor} " : ''
    def prev   = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence} " : ''
    def ctx    = context.name != 'NO_FILE' ? "--context ${context} " : ''
    """
    mhcmatch rank pairs ${table} \\
        --cls ${cls} \\
        --tier ${params.mhcmatch_tier} \\
        --passthrough --prefix '${params.mhcmatch_rerank_prefix}' \\
        --epitope ${epi} \\
        ${keepArgs()}${ctx}${tumor}${prev}${rankExtras()}${task.ext.args ?: ''} \\
        --out ${prefix}.${cls}.epitopes.mhcmatch.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The caller's columns lead and a stub cannot know them, so it types what the command ADDS.
    """
    python -c "from mhcmatch import rank; print('\\t'.join('${params.mhcmatch_rerank_prefix}' + c \\
        for c in rank.columns(\\
        extended=${isOn(params.mhcmatch_rank_extended) ? 'True' : 'False'}, \\
        annotate=${isOn(params.mhcmatch_rank_annotate) ? 'True' : 'False'}, \\
        core=${isOn(params.mhcmatch_rank_core) ? 'True' : 'False'}, \\
        cls='${cls}', species='${speciesOf()}', mode='${params.mhcmatch_rank_epitope}')))" \\
        > ${prefix}.${cls}.epitopes.mhcmatch.tsv
    """
}


process MHCMATCH_NEOAG {
    tag "${meta.id}:${cls}"
    label 'process_single'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    input:
    tuple val(meta), path(peptides), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.neoag.tsv"), emit: neoag

    script:
    def prefix = armPrefix(meta, task.ext.prefix)
    """
    mhcmatch neoag --peptides ${peptides} --cls ${cls} \\
        --max-subs ${params.mhcmatch_neoag_max_subs} \\
        ${isOn(params.mhcmatch_neoag_core) ? '--core ' : ''}${task.ext.args ?: ''} \\
        --out ${prefix}.${cls}.mhcmatch.neoag.tsv
    """

    stub:
    """
    python -c "
from mhcmatch.mimicry import NEOAG_COLUMNS
from mhcmatch.rank import CORE_COLUMNS
print('\\t'.join(('peptide',) + NEOAG_COLUMNS + (CORE_COLUMNS if ${isOn(params.mhcmatch_neoag_core) ? 'True' : 'False'} else ())))" \\
        > ${armPrefix(meta, task.ext.prefix)}.${cls}.mhcmatch.neoag.tsv
    """
}


process MHCMATCH_MIMICRY {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    input:
    tuple val(meta), path(peptides), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.mimicry.tsv"), emit: mimicry

    script:
    def prefix = armPrefix(meta, task.ext.prefix)
    """
    mhcmatch mimicry --peptides ${peptides} --cls ${cls} \\
        ${isOn(params.mhcmatch_mimicry_annotate) ? '--annotate ' : ''}${task.ext.args ?: ''} \\
        --out ${prefix}.${cls}.mhcmatch.mimicry.tsv
    """

    stub:
    """
    python -c "
from mhcmatch.rank import MIMICRY_PAIRS
print('\\t'.join(['peptide', 'logodds', 'autoimmune'] + [f'{c}_{ch}' for c, ch in MIMICRY_PAIRS]))" \\
        > ${armPrefix(meta, task.ext.prefix)}.${cls}.mhcmatch.mimicry.tsv
    """
}


process MHCMATCH_CASSETTE_SELECT {
    tag "${meta.id}"
    label 'process_medium'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    // Fixed k, where MHCMATCH_CASSETTE sizes by the per-allotype stopping rule of `--n0`. Both are
    // real answers to "how many units": `--n0` says how many the recipient's allotypes can carry,
    // `-k` says how many will be manufactured. A trial with a committed construct size needs `-k`.
    //
    // NO `--species` and no `--tier`: `cassette select` accepts neither and exits 2 if handed one.
    // The selection is over a scored pool and reads no panel.
    input:
    tuple val(meta), path(candidates), val(alleles)

    output:
    tuple val(meta), path("*.vaccine.units.tsv"), emit: units

    script:
    def prefix = armPrefix(meta, task.ext.prefix)
    def tol    = params.mhcmatch_cassette_tol ? "--tol ${params.mhcmatch_cassette_tol} " : ''
    def scol   = scoreColumn(meta)
    def prev   = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence} " : ''
    def rho    = params.mhcmatch_cassette_rho ? "--rho ${params.mhcmatch_cassette_rho} " : ''
    // The donor's DISTINCT allotypes: the denominator coverage is reported against, so an allotype
    // holding zero units is visible. Without it coverage is taken over the labels the cassette
    // happens to carry and cannot see the one it missed.
    def uni    = alleles ? "--universe '${alleles}' " : ''
    """
    mhcmatch cassette select \\
        --candidates ${candidates} \\
        -k ${params.mhcmatch_cassette_k} \\
        --block-live ${params.mhcmatch_hla_loss} \\
        ${tol}${scol ? "--score-column ${scol} " : ''}${prev}${rho}${uni}${task.ext.args ?: ''} \\
        --passthrough \\
        --out ${prefix}.vaccine.units.tsv
    """

    stub:
    // Asks the library, like every other stub here. Typed by hand it read 23 columns against a real
    // 25. `group` is dropped because a real run emits it only under `--group-column`.
    """
    python -c "from mhcmatch.cassette import SELECT_COLUMNS as C; print('\\t'.join(c for c in C if c != 'group'))" \\
        > ${armPrefix(meta, task.ext.prefix)}.vaccine.units.tsv
    """
}


process MHCMATCH_CASSETTE {
    tag "${meta.id}"
    label 'process_high'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    // `alleles` is the class-I list as a String, or a `[mhc1: '...', mhc2: '...']` Map when the
    // caller has a per-donor class-II list too. The Map is how a per-donor class-II list gets here
    // without a sixth tuple element; the class-II half is read for the cassette map's `self_help`
    // (does this unit's CD8 epitope have CD4 help from the SAME unit?) and nowhere else.
    input:
    tuple val(meta), path(candidates), path(context), val(alleles), val(cls)

    output:
    tuple val(meta), path("*.cassette.tsv"),      emit: report
    tuple val(meta), path("*.cassette.faa"),      emit: protein
    tuple val(meta), path("*.cassette.fna"),      emit: cds
    tuple val(meta), path("*.cassette.map.tsv"),  optional: true, emit: map
    tuple val(meta), path("*.cassette.map.json"), optional: true, emit: map_json

    script:
    def prefix = armPrefix(meta, task.ext.prefix)
    // **`build` selects; `order` lays out what it was handed.** Both arms put CASSETTE_SELECT in
    // front, which has already chosen exactly `-k` units -- so `build` re-selects them under the
    // per-allotype `--n0` rule and throws most away (measured on a real donor: `-k 20` in, 2 units
    // and 54 aa out). `order` is the same code path with the sizing rule skipped, so the safety
    // screen still runs and the layout, spacer and map are byte-for-byte what a selection gets.
    def verb   = task.ext.verb ?: 'order'
    // **The off state announces itself, on every task.** A missing safety check that says nothing
    // looks exactly like one that ran and found nothing.
    def screen = isOn(params.mhcmatch_cassette_screen) ? '--screen ' : ''
    def warn   = screen ? '' : ("echo '# NO SAFETY SCREEN RAN: --mhcmatch_cassette_screen is off, " +
                 "so no unit was withdrawn for essential-tissue self-origin and this cassette " +
                 "carries whatever it was handed.' >&2")
    def ctx    = context.name != 'NO_FILE' ? "--context ${context} " : ''
    // The long window WITHOUT a context FASTA: a scored candidate table often already carries it in
    // a column of its own (the fixtures spell it `context_peptide`; every upstream spells it
    // differently, which is why this is a parameter and not a default). Ignored when both are set.
    def ucol   = (context.name == 'NO_FILE' && params.mhcmatch_cassette_unit_column)
                     ? "--unit-column ${params.mhcmatch_cassette_unit_column} " : ''
    def a1     = alleles instanceof Map ? (alleles.mhc1 ?: '') : alleles
    def a2     = (alleles instanceof Map ? alleles.mhc2 : null) ?: params.mhcmatch_cassette_map_alleles_mhc2
    // **The two classes do not share a %rank cut-off**, which is why this passes a NetMHCpan TIER
    // and not a number: class I is strong <= 0.5 / weak <= 2.0, class II strong <= 2.0 / weak
    // <= 10.0. A flat `--map-threshold 2.0` is therefore the weak cut for class I and the STRONG
    // cut for class II, and it is why one mouse construct reported zero class-II epitopes with its
    // best window at %rank 4.095.
    def mt1    = params.mhcmatch_cassette_map_threshold      ? "--map-threshold ${params.mhcmatch_cassette_map_threshold} " : ''
    def mt2    = params.mhcmatch_cassette_map_threshold_mhc2 ? "--map-threshold-mhc2 ${params.mhcmatch_cassette_map_threshold_mhc2} " : ''
    def mapArg = isOn(params.mhcmatch_cassette_map)
                     ? "--map ${prefix}.cassette.map.tsv --map-json ${prefix}.cassette.map.json " +
                       "--map-binder ${params.mhcmatch_cassette_map_binder} " +
                       "${mt1}${mt2}${a2 ? "--map-alleles-mhc2 '${a2}' " : ''}" : ''
    // Quota composition: fill declared slot budgets so that at least k of each arm is expected to
    // respond, rather than taking the ranked top. Off unless a quota is given -- the arms and their
    // targets are a trial-design decision with no defensible default. With a quota the FASTAs carry
    // TWO records (`cassette_composed` and `cassette_topk`), so the comparison is on the
    // recipient's own candidates; the map describes the composed one.
    def quota  = params.mhcmatch_cassette_quota
                     ? "--quota '${params.mhcmatch_cassette_quota}' " +
                       "--block-live ${params.mhcmatch_quota_block_live} " +
                       "--evenness ${params.mhcmatch_cassette_evenness} " : ''
    // `order` does not size, so it needs no capacity estimate.
    def n0arg  = verb == 'order' ? '' : "--n0 ${params.mhcmatch_cassette_n0} "
    if (verb != 'order' && params.mhcmatch_cassette_n0 == null)
        error "params.mhcmatch_cassette_n0 is required by `cassette build` and has no default: " +
              "per-allotype capacity is not fitted by anything in the public record, so the value " +
              "is yours to set and it is recorded in the output"
    // **A unit is the long (~27 aa) window around the mutation, and there are exactly two sources
    // for it.** With neither, `_read_units` falls back to `peptide` -- the MINIMAL epitope, and a
    // 9-mer loads onto any cell without costimulation, which is the tolerising configuration
    // (PMID 17911588). So this refuses rather than defaulting.
    if (!ctx && !ucol)
        error "MHCMATCH_CASSETTE (${meta.id}): no source for the cassette unit. Give this sample a " +
              "`windows` file in the samplesheet (it becomes --context), or set " +
              "--mhcmatch_cassette_unit_column to the column of your candidate table holding the " +
              "long window (the fixtures spell it `context_peptide`). Falling back to the minimal " +
              "epitope is not a smaller version of the right thing -- it is the tolerising " +
              "configuration (PMID 17911588) -- so it is refused rather than defaulted"
    """
    ${warn}
    mhcmatch cassette ${verb} \\
        --candidates ${candidates} \\
        --alleles '${a1}' \\
        --cls ${cls} \\
        --tier ${params.mhcmatch_tier} \\
        ${ctx}${ucol}${n0arg}${screen}${quota}${mapArg}${task.ext.args ?: ''} \\
        --fasta ${prefix}.cassette.faa \\
        --fasta-nt ${prefix}.cassette.fna \\
        --out ${prefix}.cassette.tsv
    """

    stub:
    def prefix = armPrefix(meta, task.ext.prefix)
    """
    printf 'section\\ti\\tkey\\tvalue\\tdetail\\n'   > ${prefix}.cassette.tsv
    printf '>cassette units=0 spacer=null\\n\\n'     > ${prefix}.cassette.faa
    printf '>cassette_cds units=0 spacer=null\\n\\n' > ${prefix}.cassette.fna
    python -c "from mhcmatch.vector import MAP_COLUMNS as C; print('\\t'.join(C))" \\
        > ${prefix}.cassette.map.tsv
    echo '{"units": [], "features": []}' > ${prefix}.cassette.map.json
    """
}


process MHCMATCH_CASSETTE_SCORE {
    tag "${arm}:${tables.size()}"
    label 'process_single'

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere, and
    // a config file has no equivalent -- which is why this one directive is not in
    // nextflow.config beside `container`.
    conda "${moduleDir}/environment.yml"

    // `moduleDir`, not `projectDir`: an integrator's entry script is elsewhere.
    conda \"${moduleDir}/environment.yml\"

    // THE WHOLE POINT OF THIS PROCESS IS THAT IT IS NOT PER DONOR.
    //
    // `rank` anchors `p_response` on the batch it is handed, so a per-donor invocation makes every
    // donor's mean candidate probability equal the declared prevalence, whatever their pool holds.
    // Measured on 7,261 TCGA donors with pools spanning 1 to 5,221 candidates: every per-donor-
    // anchored pool mean lands on 0.060163, standard deviation 2.75e-17. Two donors' numbers are
    // then the same number and a cross-donor triage reads noise.
    //
    // `cassette score` fits ONE offset over every row it is given, so this takes the COLLECTED
    // tables of the whole run -- which is why the input is a plain `path` list with no `meta`.
    //
    // `tables` is the **units** table `cassette select` writes (`donor, slot, peptide, allele,
    // gene, score, p`), NOT the `.cassette.tsv` report: that report is long-form with the peptide
    // absent and `p` inside a free-text `detail` field, so `cassette score` could never read it.
    input:
    tuple val(arm), path(tables), path(pools)

    output:
    tuple val(arm), path("*.cassette_score.tsv"), emit: score

    script:
    def prefix = task.ext.prefix ?: (arm ? "cohort.${arm}" : 'cohort')
    def prev   = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence} " : ''
    def rho    = params.mhcmatch_cassette_rho ? "--rho ${params.mhcmatch_cassette_rho} " : ''
    def per    = isOn(params.mhcmatch_cassette_per_donor_offset) ? '--per-donor-offset ' : ''
    // `pools` is a collected List, so `.name` is a List of names and comparing it to a String is
    // always false. Ask whether any real file arrived instead, which is what was meant.
    def havePool = pools instanceof List ? pools.any { it.name != 'NO_FILE' } : pools.name != 'NO_FILE'
    // Which column carries the aggregate, resolved per ARM: on the rerank arm the POOL is the
    // caller's own table and HAS a `score` column of theirs, so an unqualified fallback scores the
    // cassette on ours and the pool on theirs -- and `lam` compares the two.
    def scol   = scoreColumn([arm: arm])
    """
    # **One call, every donor's file.** `cassette score` takes `--cassettes` and `--pool` as file
    # LISTS, reads each with its own header, and falls back to the basename up to the first dot when
    # a file's `donor` column is absent. This replaced embedded awk that reprojected every table by
    # name and was wrong twice -- once applying the first file's header to a second file's rows (two
    # mouse lines differing by two columns, giving a `gene_name` holding A/C/G/T), once re-emitting
    # the header per file so every table after the first contributed a header as a data row.
    mhcmatch cassette score \\
        --cassettes ${tables} \\
        --block-live ${params.mhcmatch_hla_loss} \\
        ${havePool ? "--pool ${pools} " : ''}${prev}${rho}${per}${scol ? "--score-column ${scol} " : ''}${task.ext.args ?: ''} \\
        --out ${prefix}.cassette_score.tsv
    """

    stub:
    """
    python -c "from mhcmatch.cassette import SCORE_COLUMNS as C; print('\\t'.join(C))" \\
        > ${task.ext.prefix ?: (arm ? "cohort.${arm}" : 'cohort')}.cassette_score.tsv
    """
}
