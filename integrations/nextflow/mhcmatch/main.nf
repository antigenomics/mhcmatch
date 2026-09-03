// mhcmatch as a set of drop-in nf-core-style local processes.
//
// Nine processes, in pipeline order:
//
//   MHCMATCH_ALLELES   an HLA typing file -> the allele list every other process takes
//   MHCMATCH_PREDICT   variant windows -> per-allele presentation, affinity, agretopicity
//   MHCMATCH_RANK      candidates      -> the fitted EPIC aggregate, one ordered table
//                                        (carries `occupancy` alongside `agretopicity`)
//   MHCMATCH_RERANK    a caller's OWN candidate table -> the same aggregate, appended to it
//   MHCMATCH_NEOAG     peptides        -> proximity to the tested-neoantigen database
//   MHCMATCH_MIMICRY   peptides        -> the six signed self/viral/thymus channels + autoimmune
//   MHCMATCH_CASSETTE_SELECT  a pool   -> the k units to manufacture (default k = 20)
//   MHCMATCH_CASSETTE  ranked units    -> a screened polyepitope cassette, amino acid and CDS
//   MHCMATCH_CASSETTE_SCORE  every donor -> one shared calibration, so donors are comparable
//
// Two entry points, and they are different objects. `subworkflows/*.nf` are for a pipeline that
// `include`s these processes into its own channel topology; `pipeline.nf` is a runnable workflow
// over a directory of files, for a caller who wants the chain and not the wiring.
//
// `subworkflows/mhcmatch.nf` chains them; see ./README.md for the input and output contract of each.
//
// Two conventions worth knowing before editing:
//
//  * **No stub types a header.** Every stub asks the installed library for its own schema, because a
//    hand-copied header is one that drifts -- this module shipped an 18-column `scored.csv` stub
//    against a 57-column real table until 2026-08-18. `predict.SCORED_COLUMNS`,
//    `predict.NATIVE_COLUMNS`, `rank.columns()` and `mimicry.NEOAG_COLUMNS` are the sources of
//    truth. What a stub cannot know is the caller's own columns: `neoag` and `mimicry` carry every
//    non-`peptide` column of a `--peptides` TSV straight through, so a real run fed `ranked.tsv`
//    emits those ahead of the ones below. A stub types the schema the command *adds*.
//  * **Species follows `params.genome`**, mapped in nextflow.config via ext.args exactly as the ARDA
//    module does, so there is no extra parameter to configure.

//: Is a boolean parameter on?
//:
//: **`--some_flag false` on the command line arrives as the STRING "false", which is truthy in
//: Groovy**, so the plain `params.x ? '--flag' : ''` idiom passes the flag a user just tried to
//: disable. And it cannot be fixed in `nextflow.config`: a config statement is evaluated when the
//: config is parsed, which is BEFORE Nextflow applies `--param` from the command line, so a
//: coercion written there is silently overwritten by the very value it exists to coerce. It has to
//: happen at the point of use -- here, and inside the resource closures in the two config files,
//: which are evaluated per task and therefore late enough.
//:
//: The direction that matters is the reverse one: somebody who believes they enabled
//: `--mhcmatch_vector_screen` and did not gets a cassette with no safety check and no error.
def isOn(v) {
    v != null && !(v.toString().toLowerCase() in ['false', '0', 'no', ''])
}


process MHCMATCH_PREDICT {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    input:
    tuple val(meta), path(fasta), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.scored.csv"), emit: scored
    // `native_tsv`, not `native`: `native` is a Java/Groovy reserved word (a method modifier), and a
    // parser that accepts it as an `emit:` name is doing us a favour rather than following a rule.
    // Nextflow 21.10.6 does not -- it fails the WHOLE module with "Unexpected input: '{'" pointing at
    // the enclosing `process {`, which reads as a corrupt file rather than one bad identifier.
    // Measured on Aldan-3 2026-08-23: `emit: native` fails as a tuple output, as a `path` output and
    // as the sole output; `emit: nativ` and `emit: native_tsv` all parse.
    tuple val(meta), val(cls), path("*.mhcmatch.native.tsv"), emit: native_tsv
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tier   = params.mhcmatch_tier ?: 'full'
    def rank   = params.mhcmatch_rank_threshold ?: 2.0
    def core   = isOn(params.mhcmatch_predict_core) ? '--core ' : ''
    """
    mhcmatch predict ${fasta} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${tier} \\
        --rank-threshold ${rank} \\
        ${core}${args} \\
        --scored-csv ${prefix}.${cls}.mhcmatch.scored.csv \\
        --native ${prefix}.${cls}.mhcmatch.native.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c "from mhcmatch.predict import SCORED_COLUMNS as C; print(','.join(C))" \\
        > ${prefix}.${cls}.mhcmatch.scored.csv
    python -c "from mhcmatch.predict import NATIVE_COLUMNS as C; print('\\t'.join(C))" \\
        > ${prefix}.${cls}.mhcmatch.native.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_RANK {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    // `rank` reads the known-epitope sets, the mimicry references and the expression tables on top
    // of the ligand panel. The image bakes them (`bootstrap --reference`); a bare `bootstrap` image
    // will reach for HuggingFace from the compute node instead.
    input:
    tuple val(meta), path(input), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.ranked.tsv"), emit: ranked
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode   = params.mhcmatch_rank_mode ?: 'fasta'
    def tier   = params.mhcmatch_tier ?: 'full'
    def tumor  = params.mhcmatch_tumor ? "--tumor ${params.mhcmatch_tumor}" : ''
    def score  = params.mhcmatch_rank_score ? "--score ${params.mhcmatch_rank_score} " : ''
    // The pool prevalence `p_response` is anchored on. It is a PRIOR about this cohort's candidate
    // list, not a model output, so it is a pipeline parameter rather than a default buried in a
    // process. Left unset, the CLI uses TESLA's 37 of 615.
    def prev   = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence} " : ''
    def extra  = (isOn(params.mhcmatch_rank_extended) ? '--extended ' : '') +
                 (isOn(params.mhcmatch_rank_annotate) ? '--annotate ' : '') +
                 (isOn(params.mhcmatch_rank_core)     ? '--core '     : '')
    """
    mhcmatch rank ${mode} ${input} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${tier} \\
        ${tumor} ${score}${prev}${extra}${args} \\
        --out ${prefix}.${cls}.mhcmatch.ranked.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext    = isOn(params.mhcmatch_rank_extended) ? 'True' : 'False'
    def ann    = isOn(params.mhcmatch_rank_annotate) ? 'True' : 'False'
    def sc     = params.mhcmatch_rank_score ?: 'aggregate'
    def cor    = isOn(params.mhcmatch_rank_core) ? 'True' : 'False'
    """
    python -c "from mhcmatch import rank; print('\\t'.join(rank.columns(extended=${ext}, annotate=${ann}, score='${sc}', core=${cor})))" \\
        > ${prefix}.${cls}.mhcmatch.ranked.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_NEOAG {
    tag "${meta.id}:${cls}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    input:
    tuple val(meta), path(peptides), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.neoag.tsv"), emit: neoag
    path "versions.yml",                                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def subs   = params.mhcmatch_neoag_max_subs ?: 2
    def core   = isOn(params.mhcmatch_neoag_core) ? '--core ' : ''
    """
    mhcmatch neoag --peptides ${peptides} --cls ${cls} --max-subs ${subs} ${core}${args} \\
        --out ${prefix}.${cls}.mhcmatch.neoag.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cor    = isOn(params.mhcmatch_neoag_core) ? 'True' : 'False'
    """
    python -c "
from mhcmatch.mimicry import NEOAG_COLUMNS
from mhcmatch.rank import CORE_COLUMNS
print('\\t'.join(('peptide',) + NEOAG_COLUMNS + (CORE_COLUMNS if ${cor} else ())))" > ${prefix}.${cls}.mhcmatch.neoag.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_MIMICRY {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    input:
    tuple val(meta), path(peptides), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.mimicry.tsv"), emit: mimicry
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ann    = isOn(params.mhcmatch_mimicry_annotate) ? '--annotate' : ''
    """
    mhcmatch mimicry --peptides ${peptides} --cls ${cls} ${ann} ${args} \\
        --out ${prefix}.${cls}.mhcmatch.mimicry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python -c "
from mhcmatch.rank import MIMICRY_PAIRS
cols = ['peptide', 'logodds', 'autoimmune'] + [f'{c}_{ch}' for c, ch in MIMICRY_PAIRS]
print('\\t'.join(cols))" > ${prefix}.${cls}.mhcmatch.mimicry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_CASSETTE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    // `--screen` builds one whole-proteome index per register length: ~12 GB peak each and a few
    // minutes apiece, which is why this process carries `process_high` and why the flag is a param.
    // WITHOUT IT NO SAFETY CHECK RUNS AT ALL and the cassette carries whatever it was handed.
    //
    // `alleles` is the class-I list as a String, or a `[mhc1: '...', mhc2: '...']` Map when the
    // caller has a per-donor class-II list too -- see the comment at `def a2` in the script block.
    input:
    tuple val(meta), path(candidates), path(context), val(alleles), val(cls)

    output:
    tuple val(meta), path("*.cassette.tsv"),    emit: report
    tuple val(meta), path("*.cassette.faa"),    emit: protein
    tuple val(meta), path("*.cassette.fna"),    emit: cds
    tuple val(meta), path("*.cassette.map.tsv"),  optional: true, emit: map
    tuple val(meta), path("*.cassette.map.json"), optional: true, emit: map_json
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // **`build` selects; `order` lays out what it was handed.** Both arms put
    // MHCMATCH_CASSETTE_SELECT in front, which has already chosen exactly `-k` units -- so running
    // `build` here re-selects them under the per-allotype `--n0` stopping rule and throws most
    // away. Measured on a real donor: `-k 20` in, `units=2` and 54 aa out. `order` is the same code
    // path with the sizing rule skipped (`cmd_vector`'s `order_only`), so the safety screen still
    // runs and the layout, spacer and map are byte-for-byte what a selected cassette gets.
    //
    // `build` stays the default, because `subworkflows/mhcmatch.nf` has no selector and the `--n0`
    // rule IS its answer to "how many units".
    def verb   = task.ext.verb ?: 'build'
    def n0     = params.mhcmatch_vector_n0
    def screen = isOn(params.mhcmatch_vector_screen) ? '--screen' : ''
    def ctx    = context.name != 'NO_FILE' ? "--context ${context}" : ''
    // The long window WITHOUT a context FASTA: a reranked candidate table already carries it in a
    // column (`epitope_context`, 27 aa, which is `--unit-length` exactly). `--context` is for the
    // de novo arm, where the table has minimal epitopes and the window has to be rebuilt; this is
    // for the rerank arm, where there is no window FASTA to rebuild from. Ignored when both are set.
    def ucol   = (context.name == 'NO_FILE' && params.mhcmatch_vector_unit_column)
                     ? "--unit-column ${params.mhcmatch_vector_unit_column}" : ''
    // The cassette MAP: one row per unit, linker and predicted epitope in 1-based coordinates.
    // The recipient's class-II allotypes are what make it worth more than a coordinate listing --
    // without them the map carries class I only and `self_help` (does this unit's CD8 epitope have
    // CD4 help from the SAME unit?) is not computed at all.
    //
    // **`val(alleles)` is a String OR a Map, and that is how a PER-DONOR class-II list gets here**
    // without a sixth tuple element -- which would be a breaking change to a signature
    // `subworkflows/mhcmatch.nf` and outside pipelines `include`. A String is the class-I list and
    // means exactly what it always meant; a `[mhc1: '...', mhc2: '...']` Map carries both, and
    // `subworkflows/rerank.nf` and `denovo.nf` build one per donor. A caller who passes the String
    // gets a byte-identical command line.
    //
    // `params.mhcmatch_vector_map_alleles_mhc2` stays, as the ONE literal panel a run may use for
    // every sample -- the mouse case, and any cohort typed once -- and answers when the tuple
    // carries no class-II list. Read here rather than in the config because a config statement runs
    // before Nextflow applies `--param`.
    def a1     = alleles instanceof Map ? (alleles.mhc1 ?: '') : alleles
    def a2     = (alleles instanceof Map ? alleles.mhc2 : null) ?:
                     params.mhcmatch_vector_map_alleles_mhc2
    def map2   = a2 ? "--map-alleles-mhc2 '${a2}'" : ''
    // **The two classes do not share a %rank cut-off**, which is why this passes a NetMHCpan
    // *tier* and not a number. NetMHCpan calls class I strong at <= 0.5 and weak at <= 2.0;
    // NetMHCIIpan calls class II strong at <= 2.0 and weak at <= 10.0. The old `--map-threshold 2.0`
    // was therefore the weak cut for class I and the STRONG cut for class II, and it is why one
    // mouse construct reported zero class-II epitopes with its best window at %rank 4.095.
    // `weak` is the default because the map reports and selects nothing.
    def mbind = params.mhcmatch_vector_map_binder ?: 'weak'
    def mt1   = params.mhcmatch_vector_map_threshold      ? "--map-threshold ${params.mhcmatch_vector_map_threshold}" : ''
    def mt2   = params.mhcmatch_vector_map_threshold_mhc2 ? "--map-threshold-mhc2 ${params.mhcmatch_vector_map_threshold_mhc2}" : ''
    def mapArg = isOn(params.mhcmatch_vector_map)
                     ? "--map ${prefix}.cassette.map.tsv --map-json ${prefix}.cassette.map.json " +
                       "--map-binder ${mbind} ${mt1} ${mt2} ${map2}" : ''
    // Quota composition: fill declared slot budgets so that at least k of each arm is expected to
    // respond, rather than taking the ranked top. Off unless a quota is given, because the arms and
    // their targets are a trial-design decision and there is no defensible default for them.
    //
    // With a quota, `.cassette.faa` / `.cassette.fna` carry TWO records -- `cassette_composed` and
    // `cassette_topk`, the same slot budgets filled by score alone -- so the comparison is on the
    // recipient's own candidates rather than asserted. The map describes the composed one. Without
    // a quota each file carries the single `cassette` record it always did.
    def quota  = params.mhcmatch_vector_quota
                     ? "--quota '${params.mhcmatch_vector_quota}' " +
                       "--block-live ${params.mhcmatch_vector_block_live ?: 0.5} " +
                       "--evenness ${params.mhcmatch_vector_evenness ?: 0.0}" : ''
    // `order` does not size, so it needs no capacity estimate.
    def n0arg  = verb == 'order' ? '' : "--n0 ${n0}"
    // The reference panel tier, which `cassette build`/`order` DO accept -- they take
    // `_add_vector_opts`, and that helper carries `--pmhc/--tier/--species` alongside the assembly
    // flags. Without this the cassette silently used `full` while `predict`, `rank` and `neoag`
    // used whatever `--mhcmatch_tier` said, so a run set to `shortlist` scored its candidates
    // against one panel and screened its cassette against another.
    def tier   = params.mhcmatch_tier ?: 'full'
    if (verb != 'order' && n0 == null) error "params.mhcmatch_vector_n0 is required and has no default: per-allotype capacity is not fitted by anything in the public record, so the value is yours to set and it is recorded in the output"
    """
    mhcmatch cassette ${verb} \\
        --candidates ${candidates} ${ctx} ${ucol} \\
        ${n0arg} \\
        --alleles '${a1}' \\
        --cls ${cls} \\
        --tier ${tier} \\
        ${screen} ${quota} ${mapArg} ${args} \\
        --fasta ${prefix}.cassette.faa \\
        --fasta-nt ${prefix}.cassette.fna \\
        --out ${prefix}.cassette.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'section\\ti\\tkey\\tvalue\\tdetail\\n' > ${prefix}.cassette.tsv
    printf '>cassette units=0 spacer=null\\n\\n'   > ${prefix}.cassette.faa
    printf '>cassette_cds units=0 spacer=null\\n\\n' > ${prefix}.cassette.fna
    python -c "from mhcmatch.vector import MAP_COLUMNS as C; print('\\t'.join(C))" \\
        > ${prefix}.cassette.map.tsv
    echo '{"units": [], "features": []}' > ${prefix}.cassette.map.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_CASSETTE_SCORE {
    tag "cohort:${tables.size()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    // THE WHOLE POINT OF THIS PROCESS IS THAT IT IS NOT PER DONOR.
    //
    // `rank` runs once per sample, and `p_response` is anchored on the batch it is handed -- so a
    // per-donor invocation makes every donor's mean candidate probability equal the declared
    // prevalence, whatever their pool holds. Measured on 7,261 TCGA donors with pools spanning 1 to
    // 5,221 candidates: every per-donor-anchored pool mean lands on 0.060163, standard deviation
    // 2.75e-17. Two donors' numbers are then the same number and a cross-donor triage reads noise.
    //
    // `cassette score` fits ONE offset over every row it is given, so this process takes the
    // COLLECTED tables of the whole run. That is why its input is a plain `path` list and not a
    // `tuple val(meta), ...`: there is no per-sample meta, because the calibration is the cohort's.
    //
    // Set `params.mhcmatch_cassette_per_donor_offset` only if you want the ENRICHMENT reading --
    // how far each donor's chosen units sit above their own background -- and know that the result
    // is no longer a probability and no longer comparable between donors.
    // `tables` is the **units** table `cassette select` writes -- `donor, slot, peptide, allele,
    // gene, score, p` -- and NOT the `.cassette.tsv` report. That report is long-form
    // (`section, i, key, value, detail`) with the peptide absent and `p` inside a free-text
    // `detail` field, so `cassette score` could never read it: it wants one row per manufactured
    // unit. Every subworkflow passed the report here until 2026-09-02, which means this process had
    // never once completed.
    input:
    path tables
    path pools

    output:
    path "*.cassette_score.tsv", emit: score
    path "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prev = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence}" : ''
    def rho  = params.mhcmatch_cassette_rho ? "--rho ${params.mhcmatch_cassette_rho}" : ''
    def per  = isOn(params.mhcmatch_cassette_per_donor_offset) ? '--per-donor-offset' : ''
    def pool = pools.name != 'NO_FILE' ? "--pool cohort.pool.tsv" : ''
    // Which column carries the aggregate. `_cassette_rows` resolves `score` / `aggregate` / `epic`
    // when not told, and on the rerank arm the POOL is the caller's own table, which has a `score`
    // column of theirs -- so the cassette would be scored on ours and the pool on theirs, and
    // `lam` compares the two. Naming the prefixed one is right for both files: the pool has it, and
    // the units table does not, so that one falls back to its own `score` and says which it used.
    // Which column holds the aggregate, resolved per ARM: on the rerank arm the pool is the
    // caller's own table and HAS a `score` column of theirs, so an unqualified fallback scores the
    // cassette on ours and the pool on theirs -- and `lam` compares the two.
    def scoreCol = task.ext.score_column ?: (params.mhcmatch_cassette_score_column ?: '')
    def scol = scoreCol ? "--score-column ${scoreCol}" : ''
    // Cohort-level and therefore no `meta`, but still one file per ARM into one publishDir.
    def prefix = task.ext.prefix ?: 'cohort'
    """
    # **Projected to four columns by NAME, not concatenated whole.** `cassette score` needs exactly
    # `donor, peptide, allele, score`, and pasting whole tables together assumes every sample has
    # the same schema -- which under `--passthrough` is false by construction, because the caller's
    # own columns travel with each sample and two samples may come from different upstream runs.
    # Measured on two mouse lines: their tables differed by two columns, so the header of the first
    # was applied to the rows of the second and every field after that point was off by two. The
    # result was not an error -- it was a `gene_name` holding `A`/`C`/`G`/`T` and a peptide column
    # that silently held something else.
    #
    # One awk over the whole file list, not one per file: a per-file loop resets awk's state, so it
    # re-emitted the header for every table after the first and those became data rows.
    project() {
        awk -F'\\t' -v OFS='\\t' -v SCOL="\$1" '
            function pick(list,   n, i, a) {
                n = split(list, a, ",")
                for (i = 1; i <= n; i++) if (a[i] != "" && (a[i] in H)) return H[a[i]]
                return 0
            }
            FNR == 1 {
                delete H; for (i = 1; i <= NF; i++) H[\$i] = i
                pc = pick("peptide,epitope")
                ac = pick("allele,mm_allele_scored,mm_allele,best_allele,allele_scored")
                sc = pick(SCOL ",score,aggregate,epic,mm_score")
                dc = pick("donor,patient")
                n = split(FILENAME, p, "/"); split(p[n], q, "."); d = q[1]
                if (!pc || !sc) {
                    printf "%s: no peptide and/or score column\\n", FILENAME > "/dev/stderr"; exit 3
                }
                if (!seen) { print "donor", "peptide", "allele", "score"; seen = 1 }
                next
            }
            {
                # The filename wins over a `donor` column holding `-`: `cassette select` writes
                # that when the pool it was given had no donor column of its own, and taking it
                # literally puts every sample in one group -- which is the cross-donor comparison
                # this process exists for.
                dv = dc ? \$dc : ""
                print (dv != "" && dv != "-" ? dv : d), \$pc, (ac ? \$ac : ""), \$sc
            }
        ' "\${@:2}"
    }
    project '${scoreCol}' ${tables} > cohort.cassettes.tsv
    if [ -n "${pool}" ]; then
        project '${scoreCol}' ${pools} > cohort.pool.tsv
    fi

    mhcmatch cassette score \\
        --cassettes cohort.cassettes.tsv \\
        ${pool} ${prev} ${rho} ${per} ${scol} ${args} \\
        --out ${prefix}.cassette_score.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: 'cohort'
    """
    printf 'donor\\tk\\tyield\\tp_mean\\tp_at_least\\toffset\\trho\\tn_effective\\tlam\\n' \\
        > ${prefix}.cassette_score.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_ALLELES {
    tag "${meta.id}:${cls}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    // **The step whose absence is silent.** Every HLA typer writes the G-group form
    // (`A*01:01:01G`) and the pseudosequence tables are keyed at two fields, so an untrimmed name
    // resolves to nothing -- and `Store._allele_set` drops what it cannot find without a word, so
    // the run scores against an EMPTY panel and exits 0. The class-II half is worse: a DP/DQ
    // molecule names both chains, so `DQA1*05:01` on its own is not a molecule at all and the two
    // rows of a typing file have to be joined.
    input:
    tuple val(meta), path(typing), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.alleles.txt"), emit: alleles
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mhcmatch alleles ${typing} --cls ${cls} ${args} --out ${prefix}.${cls}.mhcmatch.alleles.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '' > ${prefix}.${cls}.mhcmatch.alleles.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_RERANK {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    // `rank pairs --passthrough`: the caller's OWN candidate table comes back with every column it
    // arrived with, in its own order, plus this model's under `--prefix`, re-ordered by the
    // aggregate. That is not a join a caller can do afterwards -- `rank` splits a cell naming
    // several alleles and the best presenter stands for the row, so the output shares neither its
    // length nor its allele column with the input.
    //
    // `context` is the window FASTA the candidates were called on, and it is what makes
    // agretopicity and `d_occupancy` defined: a candidate table carries the mutant k-mer and
    // nothing the germline is recoverable from. Pass NO_FILE and every row is `wt_absent` --
    // correct, and a weaker model. Measured on one donor's 3,293 class-I candidates: 3,090 of the
    // 3,136 missense rows recover a wild type, every one of them differing at exactly one residue.
    input:
    tuple val(meta), path(table), path(context), val(cls)

    output:
    tuple val(meta), val(cls), path("*.epitopes.mhcmatch.tsv"), emit: reranked
    path "versions.yml",                                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tier   = params.mhcmatch_tier ?: 'full'
    def pre    = params.mhcmatch_rerank_prefix ?: 'mm_'
    def tumor  = params.mhcmatch_tumor ? "--tumor ${params.mhcmatch_tumor}" : ''
    def prev   = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence} " : ''
    def ctx    = context.name != 'NO_FILE' ? "--context ${context}" : ''
    def extra  = (isOn(params.mhcmatch_rank_extended) ? '--extended ' : '') +
                 (isOn(params.mhcmatch_rank_annotate) ? '--annotate ' : '') +
                 (isOn(params.mhcmatch_rank_core)     ? '--core '     : '')
    """
    mhcmatch rank pairs ${table} \\
        --cls ${cls} \\
        --tier ${tier} \\
        --passthrough --prefix '${pre}' \\
        ${ctx} ${tumor} ${prev}${extra}${args} \\
        --out ${prefix}.${cls}.epitopes.mhcmatch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pre    = params.mhcmatch_rerank_prefix ?: 'mm_'
    def ext    = isOn(params.mhcmatch_rank_extended) ? 'True' : 'False'
    def ann    = isOn(params.mhcmatch_rank_annotate) ? 'True' : 'False'
    def cor    = isOn(params.mhcmatch_rank_core) ? 'True' : 'False'
    """
    # The caller's columns lead and a stub cannot know them, so it types what the command ADDS --
    # asked of the library, never copied, so `-stub-run` cannot drift from the real shape.
    python -c "
from mhcmatch import rank
print('\\t'.join('${pre}' + c for c in rank.columns(extended=${ext}, annotate=${ann}, core=${cor})))" \\
        > ${prefix}.${cls}.epitopes.mhcmatch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}


process MHCMATCH_CASSETTE_SELECT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.mhcmatch_container

    // Fixed k, where MHCMATCH_CASSETTE sizes by the per-allotype stopping rule of `--n0`. Both are
    // real answers to "how many units": `--n0` says how many the recipient's allotypes can carry,
    // `-k` says how many will be manufactured. A trial that has already committed to a construct
    // size needs the second.
    //
    // NO `--species`: `cassette select` does not accept it and exits 2 if handed one, which is the
    // failure nextflow.config records for MIMICRY. The selection is over a scored pool and reads
    // no panel.
    input:
    tuple val(meta), path(candidates), val(alleles)

    output:
    tuple val(meta), path("*.vaccine.units.tsv"), emit: units
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def k      = params.mhcmatch_cassette_k ?: 20
    def tol    = params.mhcmatch_cassette_tol ? "--tol ${params.mhcmatch_cassette_tol}" : ''
    def scol   = params.mhcmatch_cassette_score_column
                     ? "--score-column ${params.mhcmatch_cassette_score_column}" : ''
    def prev   = params.mhcmatch_prevalence ? "--prevalence ${params.mhcmatch_prevalence}" : ''
    def rho    = params.mhcmatch_cassette_rho ? "--rho ${params.mhcmatch_cassette_rho}" : ''
    // **NOT `mhcmatch_vector_block_live`.** The two flags share a name and are different knobs:
    // on `cassette build --quota` it is P(a block is live) in the response model and defaults to
    // 0.5; on `cassette select` it is the HLA-LOSS rate and defaults to 1.0 (nothing is ever
    // lost). Passing 0.5 here makes any unit whose marginal p exceeds it unrepresentable and the
    // run stops -- measured on a real donor, 1 of 20 chosen units at p = 0.7782.
    def bl     = params.mhcmatch_cassette_block_live
                     ? "--block-live ${params.mhcmatch_cassette_block_live}" : ''
    // The donor's DISTINCT allotypes: the denominator coverage is reported against, so an allotype
    // holding zero units is visible. Without it coverage is taken over the labels the cassette
    // happens to carry and cannot see the one it missed.
    def uni    = alleles ? "--universe '${alleles}'" : ''
    """
    mhcmatch cassette select \\
        --candidates ${candidates} \\
        -k ${k} ${tol} ${scol} ${prev} ${rho} ${bl} ${uni} \\
        --passthrough ${args} \\
        --out ${prefix}.vaccine.units.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'donor\\tslot\\tpeptide\\tallele\\tgene\\tscore\\tp\\tk\\tpool_n\\toffset\\tenergy\\tlam\\trho\\tgamma\\tchannels\\tblock_live\\tselectivity\\trule\\tpi\\tnot_worse\\tdiversity\\tn_covered\\tn_allotypes\\n' \\
        > ${prefix}.vaccine.units.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}
