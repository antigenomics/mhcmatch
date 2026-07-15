// mhcmatch epitope predictor for the Gamaleya nextflow_vaccine pipeline.
//
// Drop-in alternative to MHCFLURRY_PREDICT_SCAN (class I) and the MHCII_BINDING subworkflow
// (class II): consumes the same (meta, .peptide.fasta, alleles) channel and emits a
// pipeline-compatible .scored.csv plus mhcmatch's richer native table.
//
// `cls` ('mhc1' | 'mhc2') is carried in the input tuple so one process serves both classes;
// instantiate it twice in the workflow (see deploy/README.md), mirroring MERGE_FASTAS_MHCI/_MHCII.
//
// Container: the deploy/Dockerfile image (mhcmatch + seqtree + baked panel; $MHCMATCH_PMHC set), so
// no --pmhc is needed. Pin it in conf/containers.config under `withName: MHCMATCH_PREDICT`.

process MHCMATCH_PREDICT {
    tag "${meta.id}:${cls}"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta), val(alleles), val(cls)

    output:
    tuple val(meta), val(cls), path("${meta.id}.${cls}.mhcmatch.scored.csv"), emit: scored
    tuple val(meta), val(cls), path("${meta.id}.${cls}.mhcmatch.native.tsv"), emit: native
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def tier = params.mhcmatch_tier ?: 'full'
    def rank = params.mhcmatch_rank_threshold ?: 2.0
    """
    mhcmatch predict ${fasta} \\
        --alleles '${alleles}' \\
        --cls ${cls} \\
        --tier ${tier} \\
        --rank-threshold ${rank} \\
        --scored-csv ${meta.id}.${cls}.mhcmatch.scored.csv \\
        --native ${meta.id}.${cls}.mhcmatch.native.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo "type,subtype,chrom,pos,gene_name,gene_id,transcript_id,uniprot_id,tpm,ffpm,epitope,epitope_context,cluster_consensus,group,best_allele,agretopicity,affinity,affinity_percentile" > ${meta.id}.${cls}.mhcmatch.scored.csv
    echo "source\tpeptide\tbest_allele\tpercent_rank\tband" > ${meta.id}.${cls}.mhcmatch.native.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcmatch: \$(python -c "import mhcmatch; print(mhcmatch.__version__)")
    END_VERSIONS
    """
}
