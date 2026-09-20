// mhcmatch as an overlay on an existing neoantigen pipeline.
//
// **It contributes ZERO new scoring processes.** Everything below is an alias of the module in
// `../mhcmatch/`, and the only thing this file adds is `MHCMATCH_PREFLIGHT`, which computes nothing
// and exists to refuse to start. That is the whole argument for the overlay: a host pipeline that
// wires mhcmatch in by hand ends up with a FORK of our processes, and a fork drifts. A hand-wired
// integration we reviewed had drifted five ways: a stub hard-typing an 18-column header against a real 57 (the
// drift our own module repaired in 2026-09-20), an `emit:` name that no longer matched, an image
// pinned to mhcmatch 1.0.1, a `--rank-threshold 2.0` on its predictor, and a class-II bridge
// calling `predict` once per allele where one call takes the list. An alias cannot drift; a sixth
// divergence will.
//
// FIVE MODES, and `off` is the default, which is what makes this safe to merge before anyone has
// decided anything:
//
//   off       preflight only. `scored` is your object verbatim and `vaccine` is empty, so both
//             host wrappings are inert.
//   annotate  rerank with --passthrough, published BESIDE your table. Your ranking is untouched,
//             and so is your construct -- this mode is a no-op on every deliverable, by design.
//   rerank    as annotate, and then `score` is rewritten from this model and handed BACK to your
//             own selector. Your selector, your linkers, your length budget, our ordering.
//   denovo    our de novo arm on your merged window FASTAs -- a third table, orthogonal to yours.
//   cassette  rerank + select + assemble + the cohort score. Our selector replaces yours.
//
// **`annotate` and `rerank` differ by exactly one process and it is the point of the mode.** A host
// selector sorts on a column it names; publishing `mm_score` beside it changes nothing at all. See
// rescore.nf.
//
// See README.md for the lines that wire it in.

include { MHCMATCH_PREFLIGHT    } from './preflight.nf'
include { MHCMATCH_RESCORE      } from './rescore.nf'
include { MHCMATCH_RERANK_ARM   } from '../mhcmatch/subworkflows/rerank.nf'
include { MHCMATCH_DENOVO_ARM   } from '../mhcmatch/subworkflows/denovo.nf'

// A FUNCTION, not a `def x = [...]` at script scope. Nextflow 26 rejects a top-level statement
// outright -- "Statements cannot be mixed with script declarations" -- and the error names the
// include site rather than this line, so it reads as the importer being broken. A function
// declaration is a declaration and is allowed.
//: The class-I list, whichever of the two shapes the allele value has -- a bare String, or
//: `[mhc1: '...', mhc2: '...']` when the caller has the donor's class-II allotypes too. The same
//: helper the arms define, repeated here because a subworkflow cannot export one.
def mhc1Of(a) {
    a instanceof Map ? (a.mhc1 ?: '') : (a ?: '')
}

def overlayModes() {
    ['off', 'annotate', 'rerank', 'denovo', 'cassette']
}

workflow MHCMATCH_OVERLAY {

    take:
    ch_candidates      // [ val(meta), path(table), val(cls) ] -- your scored candidate table
    ch_windows         // [ val(meta), path(fasta), val(cls) ] -- your mutation-window FASTA, or empty
    ch_alleles         // [ val(meta), val(alleles) ] -- class-I list, or [mhc1:'…', mhc2:'…']

    // **There is deliberately no fourth argument.** An earlier signature took the FASTA your own
    // construct step emits, purely to hand it straight back in the four modes that do not replace
    // it. That is a channel cycle the moment a host uses BOTH seams: the overlay would need its
    // selector's OUTPUT while producing its INPUT, and Nextflow cannot schedule that. `vaccine` is
    // now empty except in `cassette` mode, and the host keeps its own with the one-line ternary in
    // README.md -- the same shape a construct-variant switch already has.

    main:
    def mode = (params.mhcmatch_overlay_mode ?: 'off').toString().toLowerCase()
    if (!overlayModes().contains(mode)) {
        error "params.mhcmatch_overlay_mode = '${mode}' is not one of ${overlayModes()}"
    }
    // Forgetting the `includeConfig` line is the failure this catches. Without it every
    // `params.mhcmatch_*` is undefined, undefined is falsy, and `mhcmatch_vector_screen` therefore
    // reads as OFF -- a cassette built with no safety screen, reported as one WARN among a dozen.
    if (mode != 'off' && !params.containsKey('mhcmatch_vector_screen')) {
        error """
        The mhcmatch params are not loaded. The overlay needs BOTH lines, not just the include:

            includeConfig '<path>/integrations/nextflow/mhcmatch/nextflow.config'
            includeConfig '<path>/integrations/nextflow/overlay/overlay.config'

        Nextflow auto-loads the config beside the ENTRY script, and this is not it.
        """.stripIndent()
    }

    // **Two records out of one assembly step is two vaccines.** With a quota the cassette FASTA
    // carries `cassette_composed` AND `cassette_topk` so the comparison is on the recipient's own
    // candidates -- useful, and not something to hand a codon optimiser, which would bracket each
    // record separately and emit two constructs under one sample name.
    if (mode == 'cassette' && params.mhcmatch_vector_quota) {
        error """
        --mhcmatch_overlay_mode cassette with --mhcmatch_vector_quota set would hand your construct
        step TWO records (cassette_composed and cassette_topk) and it would optimise both.
        Pick one: drop the quota, or run the quota comparison under mode 'rerank' and read
        ${params.outdir}/mhcmatch/*.cassette.faa yourself.
        """.stripIndent()
    }

    ch_versions = Channel.empty()

    MHCMATCH_PREFLIGHT( ch_candidates.map { meta, table, cls -> table }.first(), mode )
    ch_versions = ch_versions.mix( MHCMATCH_PREFLIGHT.out.versions )

    // Everything below is additive. `vaccine` stays EMPTY unless this run is building a cassette,
    // and `scored` is the object you passed in unless this run is rewriting it -- so in `off` mode
    // the overlay contributes nothing to either seam and both host wrappings are inert.
    ch_vaccine  = Channel.empty()
    ch_scored   = ch_candidates
    ch_reranked = Channel.empty()
    ch_denovo   = Channel.empty()
    ch_units    = Channel.empty()
    ch_score    = Channel.empty()

    if (mode in ['annotate', 'rerank', 'cassette']) {
        // `--context` is not redundancy: a candidate table carries the MUTANT k-mer and nothing the
        // germline counterpart is recoverable from, so without the window FASTA every row reads
        // `wt_absent` and agretopicity is undefined -- correct, and a weaker model.
        ch_in = ch_candidates
            .map { meta, table, cls -> [ [meta, cls], meta, table, cls ] }
            .join( ch_windows.map { meta, fa, cls -> [ [meta, cls], fa ] }, remainder: true )
            // `moduleDir` and not `projectDir`, for the reason spelled out on the same sentinel in
            // ../mhcmatch/subworkflows/rerank.nf: projectDir is the ENTRY script's directory, so a
            // host including this overlay resolved it against THEIR repo root, and test/main.nf
            // resolved it against test/../mhcmatch. Neither path exists.
            .map { key, meta, table, cls, fa ->
                [ meta, table, fa ?: file("${moduleDir}/../mhcmatch/NO_FILE"), cls ] }

        MHCMATCH_RERANK_ARM( ch_in, ch_alleles )
        ch_reranked = MHCMATCH_RERANK_ARM.out.reranked
        ch_units    = MHCMATCH_RERANK_ARM.out.units
        ch_score    = MHCMATCH_RERANK_ARM.out.score
        ch_versions = ch_versions.mix( MHCMATCH_RERANK_ARM.out.versions )

        if (mode == 'rerank') {
            // **BOTH classes, which is why this reads MHCMATCH_RERANK's own output** and not the
            // arm's class-I pool: a host construct carries class-II epitopes too and its selector
            // sorts them on the same column name.
            MHCMATCH_RESCORE( ch_reranked.map { meta, cls, tsv -> [ meta, tsv, cls ] } )
            ch_scored   = MHCMATCH_RESCORE.out.csv
            ch_versions = ch_versions.mix( MHCMATCH_RESCORE.out.versions.first() )
        }

        if (mode == 'cassette') {
            // **The one place a construct is replaced, and your own step still RUNS.** The host's
            // ternary sends the optimiser this FASTA instead of its own -- it does not delete the
            // step that would have made one. Its other outputs are what an epitope-marking step
            // and a neoantigen table read, and cutting it would take a deliverable down without
            // failing anything, which is exactly what this overlay must not do.
            ch_vaccine = MHCMATCH_RERANK_ARM.out.cassette
        }
    }

    if (mode == 'denovo') {
        // **The de novo arm takes ONE channel and carries the alleles inside it**, unlike the
        // rerank arm which takes two. That asymmetry is real: `MHCMATCH_PREDICT` needs the allele
        // list in the same tuple as the FASTA, because one process instance serves both classes and
        // the class rides in the tuple with it. Join here rather than assuming they match.
        ch_dn = ch_windows
            .map { meta, fa, cls -> [ meta, fa, cls ] }
            .combine( ch_alleles, by: 0 )
            .map { meta, fa, cls, alleles -> [ meta, fa, mhc1Of(alleles), cls ] }
        MHCMATCH_DENOVO_ARM( ch_dn )
        ch_denovo   = MHCMATCH_DENOVO_ARM.out.ranked
        ch_units    = MHCMATCH_DENOVO_ARM.out.units
        ch_score    = MHCMATCH_DENOVO_ARM.out.score
        ch_versions = ch_versions.mix( MHCMATCH_DENOVO_ARM.out.versions )
    }

    emit:
    vaccine   = ch_vaccine     // the cassette, in `cassette` mode. EMPTY in every other mode --
                               // your own construct is yours and the overlay never carries it
    scored    = ch_scored      // [ meta, table, cls ] for YOUR selector. The object you passed in,
                               // verbatim, unless mode=rerank -- then the score column comes from
                               // this model and `<column>_host` carries yours
    reranked  = ch_reranked    // your table + an `mm_` block, published beside it
    denovo    = ch_denovo
    units     = ch_units
    score     = ch_score
    preflight = MHCMATCH_PREFLIGHT.out.report
    versions  = ch_versions
}
