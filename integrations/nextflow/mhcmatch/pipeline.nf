#!/usr/bin/env nextflow
//
// mhcmatch, runnable from a directory of files.
//
//   nextflow run pipeline.nf --indir <dir> --outdir results --mode both --mhcmatch_vector_n0 8
//
// This is the **easy entry point**, not the integration surface. A pipeline that wants mhcmatch as
// a component should `include` the processes in ./main.nf or the arms in ./subworkflows/ into its
// own channel topology -- that is what those files are for, and they are unchanged by this one.
// This script exists for the caller who has files on disk and wants the chain, not the wiring.
//
// TWO ARMS, and they answer different questions:
//
//   --mode rerank   your own candidate table, re-scored and re-ordered by the EPIC aggregate, with
//                   every column you sent carried through under its own name. The deliverable is
//                   YOUR table plus a `mm_` block.
//   --mode denovo   your mutation-window FASTA, with the epitope table built entirely by mhcmatch:
//                   binding called, ranked, annotated. The deliverable is OUR table.
//   --mode both     both, independently. They do not share a cassette; each arm builds its own.
//
// Both arms end in a cassette: the k units to manufacture (`-k`, default 20) as a TSV, and the
// assembled construct as amino acids and as a CDS, with the linker chosen by minimising junctional
// binding.
//
// FILE NAMING is the entire input contract, and the sample id is the filename up to its first dot:
//
//   <id>.mhcI.epitopes.scored.tsv    /  <id>.mhcII.epitopes.scored.tsv     -> --mode rerank
//   <id>.mhcI.peptide.fasta          /  <id>.mhcII.peptide.fasta           -> --mode denovo,
//                                                                    and --context for rerank
//   <id>.alleles.tsv / <id>_norma.alleles.tsv / <id>_alleles.tsv          -> the allele list
//                    (any `<id>*alleles.tsv`; the id is what remains after the suffix)
//
// A file that does not match is ignored, and the run says how many of each it found. Pass explicit
// globs (--epitopes/--windows/--typing) when your names differ.

nextflow.enable.dsl = 2

include { MHCMATCH_ALLELES        } from './main.nf'
include { MHCMATCH_RERANK_ARM     } from './subworkflows/rerank.nf'
include { MHCMATCH_DENOVO_ARM     } from './subworkflows/denovo.nf'

// **Everything below the includes is a declaration, not a statement.** Nextflow 26.x strict syntax
// rejects a bare statement at script level -- "Statements cannot be mixed with script declarations"
// -- so `params.indir = null` and `def idOf = { ... }` both fail to compile. The params this script
// owns are declared in ./nextflow.config beside every `params.mhcmatch_*`, and these two helpers are
// FUNCTIONS rather than closures assigned to a name.

//: The sample id is the filename up to its first dot, which is the same rule
//: `MHCMATCH_CASSETTE_SCORE` already uses to derive its `donor` column.
def idOf(f) {
    f.name.replaceFirst(/\..*$/, '')
}

//: `mhcI` / `mhcII` read off the filename. Both spellings, because the pipeline schema writes the
//: first and a lot of downstream tooling writes the second.
def clsOf(f) {
    f.name.contains('.mhcII.') || f.name.contains('.mhc2.') ? 'mhc2' : 'mhc1'
}

workflow {

    if( !params.indir && !params.epitopes && !params.windows )
        error "give --indir <dir>, or explicit --epitopes / --windows globs. See ./README.md"
    if( !(params.mode in ['rerank', 'denovo', 'both']) )
        error "--mode must be rerank, denovo or both (got '${params.mode}')"

    def dir     = params.indir
    def epiGlob = params.epitopes ?: (dir ? "${dir}/*.{mhcI,mhcII}.epitopes.scored.tsv" : null)
    def winGlob = params.windows  ?: (dir ? "${dir}/*.{mhcI,mhcII}.peptide.fasta"       : null)
    def typGlob = params.typing   ?: (dir ? "${dir}/*alleles.tsv"                       : null)

    def wantRerank = params.mode in ['rerank', 'both']
    def wantDenovo = params.mode in ['denovo', 'both']

    if( wantRerank && !epiGlob ) error "--mode ${params.mode} needs --epitopes or --indir"
    if( wantDenovo && !winGlob ) error "--mode ${params.mode} needs --windows or --indir"

    // The window FASTAs serve BOTH arms: the de novo input, and the rerank arm's `--context`,
    // which is the only thing that makes agretopicity and `d_occupancy` defined there.
    ch_win = winGlob ? Channel.fromPath(winGlob, checkIfExists: false)
                              .map { f -> [ idOf(f), clsOf(f), f ] }
                     : Channel.empty()
    ch_epi = (wantRerank && epiGlob) ? Channel.fromPath(epiGlob, checkIfExists: false)
                                              .map { f -> [ idOf(f), clsOf(f), f ] }
                                     : Channel.empty()

    // Every (sample, class) this run will touch, so an allele list is produced exactly once for
    // each and both arms join against the same channel.
    ch_keys = ch_epi.map { id, cls, f -> [ id, cls ] }
                    .mix( ch_win.map { id, cls, f -> [ id, cls ] } )
                    .unique()

    // ---- the allele list, per (sample, class), as a plain string ------------------------------
    //
    // `mhcmatch alleles` is not optional plumbing. Every HLA typer writes the G-group form
    // (`A*01:01:01G`), which resolves to NO pseudosequence, and `Store._allele_set` drops what it
    // cannot find **silently** -- so a run handed a raw typing file scores against an empty panel
    // and exits 0. It also performs the join a DP/DQ heterodimer needs, since `DQA1*05:01` on its
    // own is not a molecule.
    //
    // `--alleles` / `--alleles_mhc2` bypass it with a literal used for every sample. That is the
    // mouse case and it is not a shortcut: an inbred line's H-2 haplotype is a property of the
    // line, so there is no typing file to read and nothing to infer.
    if( params.alleles || params.alleles_mhc2 ) {
        ch_alleles = ch_keys.map { id, cls ->
            [ id, cls, (cls == 'mhc2' ? params.alleles_mhc2 : params.alleles) ?: '' ] }
    }
    else if( typGlob ) {
        // **The id rule has to be as wide as the glob that feeds it.** The glob is `*alleles.tsv`,
        // which also admits `<id>_alleles.tsv` and `<id>_norma.alleles.tsv`. A narrower strip left
        // the id as `D1_alleles`, that key joined no sample, the sample lost its allele list, and
        // the de novo filter below dropped it -- warning about missing alleles, which names the
        // symptom and not the cause. One permissive suffix strip covers every spelling the glob
        // can admit, so the two can no longer disagree.
        ch_typing = Channel.fromPath(typGlob, checkIfExists: false)
                           .map { f -> [ f.name.replaceFirst(/[._-]?(norma|normal)?[._-]?alleles\.tsv$/, ''), f ] }

        // A typing file that matched the glob but joins no (sample, class) key is a naming
        // mismatch, not an absence. Say so, and name the id that was derived -- silence here is
        // exactly what turns one typo into "the panel came back empty".
        // One param and an index, NOT `{ kids, tids -> }`: `map` destructures a tuple into a
        // multi-parameter closure, `subscribe` hands the item over whole, and the two-param form
        // silently binds the entire pair to the first name.
        ch_keys.map { id, cls -> id }.unique().toList()
               .combine( ch_typing.map { id, f -> id }.unique().toList() )
               .subscribe { pair ->
                   def kids = pair[0] as List
                   def tids = pair[1] as List
                   def orphan = tids.findAll { !kids.contains(it) }
                   if( orphan )
                       log.warn "typing file(s) matched but no input sample carries that id: " +
                                "${orphan.join(', ')} -- samples seen: ${kids.join(', ')}"
               }
        MHCMATCH_ALLELES( ch_keys.combine( ch_typing, by: 0 )
                                 .map { id, cls, f -> [ [id: id, cls: cls], f, cls ] } )
        ch_alleles = MHCMATCH_ALLELES.out.alleles
                        .map { meta, cls, f -> [ meta.id, cls, f.text.trim() ] }
    }
    else {
        ch_alleles = ch_keys.map { id, cls -> [ id, cls, '' ] }
    }

    // **An empty glob is the failure this whole module keeps meeting**, so it is not allowed to be
    // silent here either: a typo in --indir, or a cohort whose filenames do not follow the
    // convention, otherwise produces a run that does nothing and exits 0.
    ch_keys.count().subscribe { n ->
        if( n == 0 )
            error "no inputs matched under --indir '${params.indir}'. Expected " +
                  "<id>.mhcI|mhcII.epitopes.scored.tsv (rerank) or <id>.mhcI|mhcII.peptide.fasta " +
                  "(de novo); pass --epitopes / --windows globs if your names differ. See ./README.md"
        else
            log.info "mhcmatch: ${n} (sample, class) input(s), mode=${params.mode}"
    }

    // ================================================================ rerank
    if( wantRerank ) {
        // `remainder: true` so a table with no window FASTA still runs: it simply has no wild
        // type, which `wt_absent` carries, and which is the honest state rather than an imputed one.
        MHCMATCH_RERANK_ARM(
            ch_epi.map { id, cls, f -> [ [id: id, cls: cls], f ] }
                  .join( ch_win.map { id, cls, f -> [ [id: id, cls: cls], f ] }, remainder: true )
                  .filter { meta, tsv, fa -> tsv != null }
                  .map { meta, tsv, fa ->
                      [ meta, tsv, fa ?: file("${moduleDir}/NO_FILE"), meta.cls ] },
            ch_alleles.filter { id, cls, a -> cls == 'mhc1' }
                      .map { id, cls, a -> [ [id: id, cls: cls], a ] }
        )
    }

    // ================================================================ de novo
    if( wantDenovo ) {
        MHCMATCH_DENOVO_ARM(
            ch_win.combine( ch_alleles, by: [0, 1] )
                  .map { id, cls, f, a -> [ [id: id, cls: cls], f, a, cls ] }
                  // No allele list, nothing to score against. `predict` REQUIRES --alleles, and a
                  // silently empty panel is the failure `mhcmatch alleles` exists to prevent -- so
                  // a sample without one is dropped loudly here rather than scored against nothing.
                  .filter { meta, f, a, cls ->
                      if( !a ) log.warn "no alleles for ${meta.id} ${cls}: skipping the de novo arm"
                      a as Boolean }
        )
    }
}
