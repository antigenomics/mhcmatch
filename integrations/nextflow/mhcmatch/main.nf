// mhcmatch as a set of drop-in nf-core-style local processes.
//
// Five processes, in pipeline order:
//
//   MHCMATCH_PREDICT   variant windows -> per-allele presentation, affinity, agretopicity
//   MHCMATCH_RANK      candidates      -> the fitted aggregate, one ordered table
//   MHCMATCH_NEOAG     peptides        -> proximity to the tested-neoantigen database
//   MHCMATCH_MIMICRY   peptides        -> the six signed self/viral/thymus channels + autoimmune
//   MHCMATCH_VECTOR    ranked units    -> a screened polyepitope cassette, amino acid and CDS
//
// `subworkflows/mhcmatch.nf` chains them; see ./README.md for the input and output contract of each.
//
// Two conventions worth knowing before editing:
//
//  * **No stub types a header.** Every stub asks the installed library for its own schema, because a
//    hand-copied header is one that drifts -- this module shipped an 18-column `scored.csv` stub
//    against a 57-column real table until 2026-08-18. `predict.SCORED_COLUMNS`,
//    `predict.NATIVE_COLUMNS` and `rank.columns()` are the sources of truth.
//  * **Species follows `params.genome`**, mapped in nextflow.config via ext.args exactly as the ARDA
//    module does, so there is no extra parameter to configure.

process MHCMATCH_PREDICT {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "mhcmatch:0.17.0"

    input:
    tuple val(meta), path(fasta), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.scored.csv"), emit: scored
    tuple val(meta), val(cls), path("*.mhcmatch.native.tsv"), emit: native
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tier   = params.mhcmatch_tier ?: 'full'
    def rank   = params.mhcmatch_rank_threshold ?: 2.0
    """
    mhcmatch predict ${fasta} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${tier} \\
        --rank-threshold ${rank} \\
        ${args} \\
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
    container "mhcmatch:0.17.0"

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
    def extra  = (params.mhcmatch_rank_extended ? '--extended ' : '') +
                 (params.mhcmatch_rank_annotate ? '--annotate ' : '')
    """
    mhcmatch rank ${mode} ${input} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${tier} \\
        ${tumor} ${extra}${args} \\
        --out ${prefix}.${cls}.mhcmatch.ranked.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext    = params.mhcmatch_rank_extended ? 'True' : 'False'
    def ann    = params.mhcmatch_rank_annotate ? 'True' : 'False'
    """
    python -c "from mhcmatch import rank; print('\\t'.join(rank.columns(extended=${ext}, annotate=${ann})))" \\
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
    container "mhcmatch:0.17.0"

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
    """
    mhcmatch neoag --peptides ${peptides} --cls ${cls} --max-subs ${subs} ${args} \\
        --out ${prefix}.${cls}.mhcmatch.neoag.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'peptide\\tneoag_distance\\tneoag_nearest\\tneoag_n_within\\tknown\\n' \\
        > ${prefix}.${cls}.mhcmatch.neoag.tsv

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
    container "mhcmatch:0.17.0"

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
    def ann    = params.mhcmatch_mimicry_annotate ? '--annotate' : ''
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


process MHCMATCH_VECTOR {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "mhcmatch:0.17.0"

    // `--screen` builds one whole-proteome index per register length: ~12 GB peak each and a few
    // minutes apiece, which is why this process carries `process_high` and why the flag is a param.
    // WITHOUT IT NO SAFETY CHECK RUNS AT ALL and the cassette carries whatever it was handed.
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
    def n0     = params.mhcmatch_vector_n0
    def screen = params.mhcmatch_vector_screen ? '--screen' : ''
    def ctx    = context.name != 'NO_FILE' ? "--context ${context}" : ''
    // The cassette MAP: one row per unit, linker and predicted epitope in 1-based coordinates.
    // `mhcmatch_vector_map_alleles_mhc2` is what makes it worth more than a coordinate listing --
    // without the recipient's class-II allotypes the map carries class I only and `self_help`
    // (does this unit's CD8 epitope have CD4 help from the SAME unit?) cannot be computed at all.
    def map2   = params.mhcmatch_vector_map_alleles_mhc2
                     ? "--map-alleles-mhc2 '${params.mhcmatch_vector_map_alleles_mhc2}'" : ''
    def mapArg = params.mhcmatch_vector_map
                     ? "--map ${prefix}.cassette.map.tsv --map-json ${prefix}.cassette.map.json " +
                       "--map-threshold ${params.mhcmatch_vector_map_threshold ?: 2.0} ${map2}" : ''
    if (n0 == null) error "params.mhcmatch_vector_n0 is required and has no default: per-allotype capacity is not fitted by anything in the public record, so the value is yours to set and it is recorded in the output"
    """
    mhcmatch vector \\
        --candidates ${candidates} ${ctx} \\
        --n0 ${n0} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        ${screen} ${mapArg} ${args} \\
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
