// mhcmatch epitope predictor, as a drop-in nf-core-style local module.
// One `mhcmatch predict` call per (sample, class); publishes to ${params.outdir}/mhcmatch/.
// See ./README.md for how to wire this into the neoantigen workflow.
//
// Drop-in for MHCFLURRY_PREDICT_SCAN (class I) and the MHCII_BINDING subworkflow (class II): consumes
// the same (meta, .peptide.fasta, alleles) channel and emits a pipeline-compatible .scored.csv plus
// mhcmatch's richer native table. `cls` ('mhc1' | 'mhc2') rides in the input tuple so one process
// serves both classes -- instantiate it twice (see ./README.md), mirroring MERGE_FASTAS_MHCI/_MHCII.
//
// Species follows the assembly the rest of the pipeline runs on: nextflow.config maps params.genome
// -> --species via ext.args (GRCh38 -> human, GRCm39 -> mouse), exactly as the ARDA module does, so
// there is no extra parameter to configure.

process MHCMATCH_PREDICT {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    // mhcmatch is pip-installable (PyPI: mhcmatch) and pulls its seqtree C++ core as a dependency.
    //   -profile conda   -> works out of the box from environment.yml
    //   -profile docker  -> build the image from the Dockerfile beside this module and push it to
    //                       your registry, then point `container` at it (or override in a config).
    conda "${moduleDir}/environment.yml"
    container "mhcmatch:0.8.0"

    input:
    tuple val(meta), path(fasta), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("*.mhcmatch.scored.csv"), emit: scored
    tuple val(meta), val(cls), path("*.mhcmatch.native.tsv"), emit: native
    path "versions.yml",                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''                 // --species (from params.genome; see nextflow.config)
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
    echo "type,subtype,chrom,pos,gene_name,gene_id,transcript_id,uniprot_id,tpm,ffpm,epitope,epitope_context,cluster_consensus,group,best_allele,agretopicity,affinity,affinity_percentile" > ${prefix}.${cls}.mhcmatch.scored.csv
    printf 'source\\tpeptide\\tbest_allele\\tpercent_rank\\tband\\n' > ${prefix}.${cls}.mhcmatch.native.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}
