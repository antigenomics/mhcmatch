// mhcmatch as an overlay on a neoantigen pipeline you already run.
//
// **It contributes ZERO new scoring processes.** Everything is the module in `../mhcmatch/`; the
// only things this directory adds are MHCMATCH_PREFLIGHT, which computes nothing and exists to
// refuse to start, and MHCMATCH_RESCORE, which writes the one column a host selector reads.
//
// That is the whole argument for an overlay: a host pipeline that wires mhcmatch in by hand ends up
// with a FORK of the processes, and a fork drifts. One we reviewed had drifted five ways within a
// few releases -- a stub hard-typing an 18-column header against a real 57, an `emit:` name that no
// longer matched, an image pinned to a stale mhcmatch, a `--rank-threshold 2.0` on its predictor (the
// weak class-I cut and the STRONG class-II one), and a class-II bridge calling `predict` once per
// allele where one call takes the list. An alias cannot drift; a sixth divergence will.
//
// FIVE MODES, and `off` is the default, which is what makes this safe to merge before anyone has
// decided anything:
//
//   off       preflight only. `scored` is your object verbatim and `vaccine` is empty, so both
//             host wrappings are inert.
//   annotate  rerank with --passthrough, published BESIDE your table. Your ranking is untouched,
//             and so is your construct -- a no-op on every deliverable, by design.
//   rerank    as annotate, and then `score` is rewritten from this model and handed BACK to your
//             own selector. Your selector, your linkers, your length budget, our ordering.
//   denovo    our de novo arm on your window FASTAs -- a third table, orthogonal to yours.
//   cassette  rerank + select + assemble + the cohort score. Our selector replaces yours.
//
// **`annotate` and `rerank` differ by exactly one process and it is the point of the mode.** A host
// selector sorts on a column it names; publishing `mm_score` beside it changes nothing at all. See
// rescore.nf.
//
// See README.md for the two lines that wire it in.

include { MHCMATCH_PREFLIGHT } from './preflight.nf'
include { MHCMATCH_RESCORE   } from './rescore.nf'
include { MHCMATCH           } from '../mhcmatch/subworkflows/mhcmatch.nf'

// A FUNCTION, not a `def x = [...]` at script scope. Nextflow 26 rejects a top-level statement
// outright -- "Statements cannot be mixed with script declarations" -- and the error names the
// INCLUDE site rather than this line, so it reads as the importer being broken.
def overlayModes() {
    ['off', 'annotate', 'rerank', 'denovo', 'cassette']
}

//: The class-I list, whichever of the two shapes the allele value has.
def mhc1Of(a) {
    a instanceof Map ? (a.mhc1 ?: '') : (a ?: '')
}

workflow MHCMATCH_OVERLAY {

    take:
    ch_candidates      // [ val(meta), path(table), val(cls) ] -- your scored candidate table
    ch_windows         // [ val(meta), path(fasta), val(cls) ] -- your window FASTA, or empty
    ch_alleles         // [ val(meta), val(alleles) ] -- class-I list, or [mhc1:'…', mhc2:'…']

    // **There is deliberately no fourth argument.** An earlier signature took the FASTA your own
    // construct step emits, purely to hand it straight back in the four modes that do not replace
    // it. That is a channel cycle the moment a host uses BOTH seams: the overlay would need its
    // selector's OUTPUT while producing its INPUT, and Nextflow cannot schedule that. `vaccine` is
    // empty except in `cassette` mode, and the host keeps its own with the one-line ternary in
    // README.md -- the same shape a construct-variant switch already has.

    main:
    def mode = (params.mhcmatch_overlay_mode ?: 'off').toString().toLowerCase()
    if( !overlayModes().contains(mode) )
        error "params.mhcmatch_overlay_mode = '${mode}' is not one of ${overlayModes()}"

    // **Two records out of one assembly step is two vaccines.** With a quota the cassette FASTA
    // carries `cassette_composed` AND `cassette_topk` so the comparison is on the recipient's own
    // candidates -- useful, and not something to hand a codon optimiser, which would bracket each
    // record separately and emit two constructs under one sample name.
    if( mode == 'cassette' && params.mhcmatch_cassette_quota )
        error "--mhcmatch_overlay_mode cassette with --mhcmatch_cassette_quota set would hand your " +
              "construct step TWO records (cassette_composed and cassette_topk) and it would " +
              "optimise both. Drop the quota, or run the quota comparison under mode 'rerank' and " +
              "read ${params.outdir}/mhcmatch/*.cassette.faa yourself."

    MHCMATCH_PREFLIGHT( ch_candidates.map { meta, table, cls -> table }.first(), mode )

    def wantRerank = mode in ['annotate', 'rerank', 'cassette']
    def arm        = wantRerank ? 'rerank' : 'denovo'

    // `--context` is not redundancy: a candidate table carries the MUTANT k-mer and nothing the
    // germline counterpart is recoverable from, so without the window FASTA every row reads
    // `wt_absent` and agretopicity is undefined -- correct, and a weaker model.
    //
    // `moduleDir` and not `projectDir` on the sentinel: projectDir is the ENTRY script's directory,
    // so a host including this overlay resolves it against THEIR repo root, where it does not exist.
    ch_rr = !wantRerank ? Channel.empty() :
        ch_candidates.map { meta, table, cls -> [ [meta, cls], meta, table, cls ] }
            .join( ch_windows.map { meta, fa, cls -> [ [meta, cls], fa ] }, remainder: true )
            .filter { key, meta, table, cls, fa -> table != null }
            .map { key, meta, table, cls, fa ->
                [ meta + [arm: 'rerank', cls: cls], table,
                  fa ?: file("${moduleDir}/../mhcmatch/NO_FILE"), cls ] }

    // The de novo arm carries the alleles INSIDE its tuple, unlike the rerank arm: MHCMATCH_PREDICT
    // needs the list beside the FASTA, because one process instance serves both classes and the
    // class rides in the tuple with it.
    ch_dn = mode != 'denovo' ? Channel.empty() :
        ch_windows.combine( ch_alleles, by: 0 )
                  .map { meta, fa, cls, alleles ->
                      [ meta + [arm: 'denovo', cls: cls], fa, mhc1Of(alleles), cls ] }

    // Keyed on the class-I meta of the arm in play, which is what the module's joins expect.
    ch_al = ch_alleles.map { meta, alleles -> [ meta + [arm: arm, cls: 'mhc1'], alleles ] }

    MHCMATCH( ch_rr, ch_dn, ch_al )

    // Everything is additive. `vaccine` stays EMPTY unless this run is building a cassette, and
    // `scored` is the object you passed in unless this run is rewriting it -- so in `off` mode the
    // overlay contributes nothing to either seam and both host wrappings are inert.
    //
    // **`rerank` mode reads MHCMATCH_RERANK's own output and not the arm's class-I pool**, because
    // a host construct carries class-II epitopes too and its selector sorts them on the same column.
    MHCMATCH_RESCORE( mode != 'rerank' ? Channel.empty()
                      : MHCMATCH.out.reranked.map { meta, cls, tsv -> [ meta, tsv, cls ] } )

    emit:
    // The cassette, in `cassette` mode; EMPTY in every other. **Your own construct step still
    // RUNS** -- the host's ternary sends the optimiser this FASTA instead of its own, it does not
    // delete the step that would have made one. That step's other outputs are what an
    // epitope-marking step and a neoantigen table read, and cutting it would take a deliverable
    // down without failing anything, which is exactly what this overlay must not do.
    vaccine   = mode == 'cassette' ? MHCMATCH.out.cassette : Channel.empty()
    // [ meta, table, cls ] for YOUR selector: the object you passed in, verbatim, unless
    // mode=rerank -- then the score column comes from this model and `<column>_host` carries yours.
    scored    = mode == 'rerank' ? MHCMATCH_RESCORE.out.csv : ch_candidates
    reranked  = MHCMATCH.out.reranked      // your table + an `mm_` block, published beside it
    denovo    = MHCMATCH.out.ranked
    units     = MHCMATCH.out.units
    score     = MHCMATCH.out.score
    preflight = MHCMATCH_PREFLIGHT.out.report
}
