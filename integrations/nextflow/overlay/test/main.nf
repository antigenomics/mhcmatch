// Stub harness for the overlay: five modes, synthetic channels shaped like a host pipeline's seam.
//
//     nextflow run test/main.nf -stub-run -c ../overlay.config --mhcmatch_overlay_mode annotate
//
// **A stub runs no command, so this cannot catch an unrecognised CLI flag** -- which is exactly the
// class that has bitten this integration before. It proves topology and nothing else. The run on a
// real fixture is the one that matters; see README.md.
nextflow.enable.dsl = 2

include { MHCMATCH_OVERLAY } from '../overlay.nf'

workflow {
    // **A bad `--table` path is invisible to a stub.** Nextflow stages a `path` input by symlink and
    // does not check the target, so a wrong relative path yields a DANGLING link in the work dir. The
    // stub never reads it and reports success; the real run then dies inside the first process with a
    // `FileNotFoundError` on a bare basename. Cost one cluster run. Check here, where it is one line.
    // The fixtures are at integrations/fixtures, which from this module is ../../fixtures.
    [ table: params.table, fasta: params.fasta ].each { name, path ->
        if (!path) error "--${name} is required: a path relative to the launch directory, or absolute"
        if (!file(path).exists()) error "--${name} ${path} does not exist (launch dir: ${launchDir})"
    }

    def meta = [ id: 'SAMPLE1' ]
    ch_cand    = Channel.of( [ meta, file(params.table),  'mhc1' ] )
    ch_windows = Channel.of( [ meta, file(params.fasta),  'mhc1' ] )
    ch_alleles = Channel.of( [ meta, 'HLA-A*02:01,HLA-B*07:02,HLA-C*07:01' ] )

    MHCMATCH_OVERLAY( ch_cand, ch_windows, ch_alleles )

    // `vaccine` carries something in exactly one mode, and nothing in the other four.
    MHCMATCH_OVERLAY.out.vaccine.view { m, f ->
        "vaccine channel carries: ${f.name}  (mode=${params.mhcmatch_overlay_mode})" }
    // `scored` is the object that went in, in every mode but `rerank`.
    MHCMATCH_OVERLAY.out.scored.view { m, f, cls ->
        "scored channel carries:  ${f.name}  (mode=${params.mhcmatch_overlay_mode})" }
}
