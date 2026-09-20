// Stub harness for the overlay: synthetic channels shaped like a host pipeline's seam.
//
//     nextflow run test/main.nf -stub-run -c ../overlay.config \
//         --mhcmatch_overlay_mode annotate \
//         --table ../../fixtures/S1.mhc1.candidates.tsv \
//         --fasta ../../fixtures/S1.mhc1.windows.fasta
//
// **A stub runs no command, so this cannot catch an unrecognised CLI flag** -- exactly the class
// that has bitten this integration before. It proves topology and nothing else; the run on a real
// fixture is the one that matters. See README.md.
nextflow.enable.dsl = 2

include { MHCMATCH_OVERLAY } from '../overlay.nf'

workflow {
    // **A bad `--table` path is invisible to a stub.** Nextflow stages a `path` input by symlink
    // and does not check the target, so a wrong relative path yields a DANGLING link in the work
    // dir: the stub never reads it and reports success, and the real run then dies inside the first
    // process with a `FileNotFoundError` on a bare basename. Cost one cluster run.
    [ table: params.table, fasta: params.fasta ].each { name, path ->
        if( !path ) error "--${name} is required: a path relative to the launch directory, or absolute"
        if( !file(path).exists() ) error "--${name} ${path} does not exist (launch dir: ${launchDir})"
    }

    def meta = [ id: 'SAMPLE1' ]
    MHCMATCH_OVERLAY(
        Channel.of( [ meta, file(params.table), 'mhc1' ] ),
        Channel.of( [ meta, file(params.fasta), 'mhc1' ] ),
        Channel.of( [ meta, 'HLA-A*02:01,HLA-B*07:02,HLA-C*07:01' ] )
    )

    // `vaccine` carries something in exactly one mode, and nothing in the other four.
    MHCMATCH_OVERLAY.out.vaccine.view { m, f ->
        "vaccine channel carries: ${f.name}  (mode=${params.mhcmatch_overlay_mode})" }
    // `scored` is the object that went in, in every mode but `rerank`.
    MHCMATCH_OVERLAY.out.scored.view { m, f, cls ->
        "scored channel carries:  ${f.name}  (mode=${params.mhcmatch_overlay_mode})" }
}
