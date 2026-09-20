#!/usr/bin/env nextflow
//
// mhcmatch, runnable from a samplesheet.
//
//   nextflow run pipeline.nf --input samplesheet.csv --outdir results --mode both
//
// This is the easy entry point, not the integration surface. A pipeline that wants mhcmatch as a
// component should `include` the processes in ./main.nf or the workflow in ./subworkflows/ into its
// own channel topology; `../overlay/` does that for a host pipeline that already reaches a
// candidate table. This script is for the caller who has files on disk and wants the chain.
//
// TWO ARMS, answering different questions:
//
//   --mode rerank   your own candidate table, re-scored and re-ordered by the EPIC aggregate, with
//                   every column you sent carried through. The deliverable is YOUR table plus an
//                   `mm_` block.
//   --mode denovo   your mutation-window FASTA, with the epitope table built entirely by mhcmatch.
//                   The deliverable is OUR table.
//   --mode both     both, independently. They do not share a cassette or a cohort calibration.
//
// THE INPUT CONTRACT IS ONE CSV, and nothing is inferred from a filename:
//
//   sample,class,candidates,windows,hla
//   S1,mhc1,S1.mhc1.candidates.tsv,S1.mhc1.windows.fasta,S1.hla.tsv
//   S1,mhc2,S1.mhc2.candidates.tsv,S1.mhc2.windows.fasta,S1.hla.tsv
//
//   sample      the id. Every output is named after it.
//   class       `mhc1` or `mhc2`, and it is not read off the file.
//   candidates  OPTIONAL -- the rerank arm's input. Any table with a peptide column and an allele
//               column; pVACseq's `*.filtered.tsv` drops in as written.
//   windows     OPTIONAL -- the de novo arm's input, AND the rerank arm's `--context`, which is
//               what makes agretopicity and `d_occupancy` defined there.
//   hla         OPTIONAL -- a typing file (OptiType wide `*_result.tsv`, HLA-LA, arcasHLA
//               `.genotype.json`, or one allele per line). Omitted for every row, `--alleles` /
//               `--alleles_mhc2` must give a literal list instead.
//
// **A relative cell resolves against the SAMPLESHEET's own directory, not the launch directory.** A
// sheet is written beside the files it names and then run from wherever the work happens; the launch
// directory is not a property of the data.
//
// The assumed upstream is the community-standard stack: nf-core/sarek -> VEP -> pVACtools for the
// variants and the candidate table, OptiType (class I) and arcasHLA / HLA-LA (class II) for the
// typing. `pvacseq generate_protein_fasta` produces the peptide-window FASTA.

nextflow.enable.dsl = 2

include { MHCMATCH_ALLELES } from './main.nf'
include { MHCMATCH         } from './subworkflows/mhcmatch.nf'

// **Everything below the includes is a declaration, not a statement.** Nextflow 26.x strict syntax
// rejects a bare statement at script level -- "Statements cannot be mixed with script declarations"
// -- so `params.input = null` and `def sheetPath = { ... }` both fail to compile. The params this
// script owns are declared in ./nextflow.config; this helper is a FUNCTION, not a closure bound to
// a name.

//: One samplesheet cell -> a staged path, or null when the cell is empty.
//:
//: `checkIfExists` is on deliberately: Nextflow stages a `path` input by symlink and does not check
//: the target, so a mistyped cell otherwise becomes a DANGLING link that `-stub-run` never reads and
//: reports success on, and the real run then dies inside a task on a bare basename.
def sheetPath(base, cell) {
    def v = cell?.toString()?.trim()
    if( !v ) return null
    file(v.startsWith('/') ? v : "${base}/${v}", checkIfExists: true)
}

workflow {

    if( !params.input )
        error "give --input <samplesheet.csv>. Columns: sample,class,candidates,windows,hla -- see ./README.md"
    if( !(params.mode in ['rerank', 'denovo', 'both']) )
        error "--mode must be rerank, denovo or both (got '${params.mode}')"

    def sheet = file(params.input, checkIfExists: true)
    def base  = sheet.parent
    def arms  = params.mode == 'both' ? ['rerank', 'denovo'] : [params.mode]

    // ---- the samplesheet, validated per row ---------------------------------------------------
    ch_rows = Channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            // **`row['class']` rather than `row.class`**, for the one column whose name is also a
            // Groovy property. Both spellings work today; the subscript form is used because it
            // cannot be read, by a reviewer or by a stricter parser, as the `getClass()` it looks
            // like, and it is what nf-core samplesheets use for the same collision.
            def id  = row.sample?.toString()?.trim()
            def cls = row['class']?.toString()?.trim()
            if( !id )
                error "samplesheet ${sheet.name}: a row has no `sample`. The columns are sample,class,candidates,windows,hla"
            if( !(cls in ['mhc1', 'mhc2']) )
                error "samplesheet ${sheet.name}, sample '${id}': `class` must be mhc1 or mhc2 (got '${cls ?: ''}')"
            def cand = sheetPath(base, row.candidates)
            def win  = sheetPath(base, row.windows)
            // A row naming neither input has nothing for either arm to read. Dropping it quietly is
            // how a cohort silently loses a donor; the sheet is the one artifact that can say so by
            // name, before any work starts.
            if( !cand && !win )
                error "samplesheet ${sheet.name}, sample '${id}' ${cls}: neither `candidates` nor " +
                      "`windows` is given, so there is nothing to score. `candidates` feeds the " +
                      "rerank arm, `windows` feeds the de novo arm and the rerank arm's --context"
            [ id, cls, cand, win, sheetPath(base, row.hla) ]
        }

    // **An empty input is the failure this module keeps meeting**, and it is not allowed to be
    // silent: a header-only sheet otherwise produces a run that does nothing and exits 0.
    ch_rows.count().subscribe { n ->
        if( n == 0 )
            error "samplesheet ${params.input} has a header and no data rows"
        else
            log.info "mhcmatch: ${n} (sample, class) row(s) from ${sheet.name}, mode=${params.mode}"
    }

    // One row per (sample, class): a duplicate row runs the whole chain twice under ONE output
    // name, which is two tasks publishing the same filename.
    //
    // One param and an index, NOT `{ key, hits -> }`: `map` destructures a tuple into a
    // multi-parameter closure, `subscribe` hands the item over WHOLE, and the two-param form
    // silently binds the entire pair to the first name.
    ch_rows.map { id, cls, cand, win, hla -> [ [id, cls], id ] }
           .groupTuple()
           .subscribe { hit ->
               if( (hit[1] as List).size() > 1 )
                   error "samplesheet ${sheet.name} names ${hit[0][0]} ${hit[0][1]} " +
                         "${(hit[1] as List).size()} times; one row per (sample, class)" }

    ch_epi  = ch_rows.filter { id, cls, cand, win, hla -> cand != null }
                     .map    { id, cls, cand, win, hla -> [ id, cls, cand ] }
    ch_win  = ch_rows.filter { id, cls, cand, win, hla -> win != null }
                     .map    { id, cls, cand, win, hla -> [ id, cls, win ] }
    ch_keys = ch_rows.map { id, cls, cand, win, hla -> [ id, cls ] }.unique()

    // A mode with no usable row is a typo in the sheet or in `--mode`, not an empty result.
    if( 'rerank' in arms )
        ch_epi.count().subscribe { n -> if( n == 0 )
            error "--mode ${params.mode} needs a `candidates` table and no row of ${sheet.name} names one" }
    if( 'denovo' in arms )
        ch_win.count().subscribe { n -> if( n == 0 )
            error "--mode ${params.mode} needs a `windows` FASTA and no row of ${sheet.name} names one" }

    // A typing file is a property of the DONOR, not of the class, so both of a sample's rows name
    // the same one and it is resolved once per sample.
    ch_typing = ch_rows.filter { id, cls, cand, win, hla -> hla != null }
                       .map    { id, cls, cand, win, hla -> [ id, hla ] }
                       .groupTuple()
                       .map { id, files ->
                           def distinct = (files as Set).toList()
                           if( distinct.size() > 1 )
                               error "samplesheet ${sheet.name}, sample '${id}': rows name " +
                                     "${distinct.size()} different `hla` files " +
                                     "(${distinct*.name.join(', ')}); a typing file is a property " +
                                     "of the donor, so every row of one sample must name the same one"
                           [ id, distinct[0] ] }

    // ---- the allele list, per (sample, class), as a plain string ------------------------------
    //
    // `mhcmatch alleles` is not optional plumbing: every HLA typer writes the G-group form
    // (`A*01:01:01G`), which resolves to NO pseudosequence, and `Store._allele_set` drops what it
    // cannot find SILENTLY -- so a run handed a raw typing file scores against an empty panel and
    // exits 0. It also performs the join a DP/DQ heterodimer needs, since `DQA1*05:01` alone is not
    // a molecule.
    //
    // `--alleles` / `--alleles_mhc2` bypass it with a literal used for every sample. That is the
    // mouse case and it is not a shortcut: an inbred line's H-2 haplotype is a property of the line,
    // so there is no typing file to read and nothing to infer.
    if( params.alleles || params.alleles_mhc2 ) {
        ch_lists = ch_keys.map { id, cls ->
            [ id, cls, (cls == 'mhc2' ? params.alleles_mhc2 : params.alleles) ?: '' ] }
    }
    else {
        ch_typing.count().subscribe { n -> if( n == 0 )
            error "no row of ${sheet.name} names an `hla` file, so there is nothing to resolve a " +
                  "panel from. Fill the `hla` column, or pass --alleles / --alleles_mhc2 with a " +
                  "literal list for every sample (the mouse case: an inbred line's haplotype is a " +
                  "property of the line)" }

        MHCMATCH_ALLELES( ch_keys.combine( ch_typing, by: 0 )
                                 .map { id, cls, f -> [ [id: id, cls: cls], f, cls ] } )
        // **The rows WITHOUT a typing file are mixed back in, carrying an empty list.** The combine
        // above is an inner join, so on its own it drops those keys entirely -- and a key missing
        // here vanishes from the de novo arm without a word.
        //
        // **`groupTuple` then prefer the non-empty, because a key can arrive twice.** `ch_typing`
        // is keyed on `id` alone, so a sample whose mhc1 row names an `hla` file and whose mhc2 row
        // leaves it blank gets a resolved list from the combine AND an empty one from the mix.
        // Without the group, whichever buffered first won -- and the empty one is emitted
        // immediately while the resolved one waits on a process -- so `self_help` was computed
        // against an empty class-II panel with nothing said.
        ch_lists = MHCMATCH_ALLELES.out.alleles
                      .map { meta, cls, f -> [ meta.id, cls, f.text.trim() ] }
                      .mix( ch_rows.filter { id, cls, cand, win, hla -> hla == null }
                                   .map    { id, cls, cand, win, hla -> [ id, cls, '' ] } )
                      .groupTuple( by: [0, 1] )
                      .map { id, cls, lists -> [ id, cls, lists.find { it } ?: '' ] }
    }

    // The donor's class-I list AND their class-II list as one value, for every arm in play.
    // MHCMATCH_CASSETTE takes `[mhc1:, mhc2:]` and computes `self_help` from the second -- whether
    // a unit's CD8 epitope has CD4 help from the SAME unit. `remainder: true` keeps a donor with no
    // class-II input, who simply gets no `self_help`; the filter drops the mirror case, an id with
    // class II and no class I, which has nothing to build a cassette from.
    ch_alleles = ch_lists.filter { id, cls, a -> cls == 'mhc1' }.map { id, cls, a -> [ id, a ] }
        .join( ch_lists.filter { id, cls, a -> cls == 'mhc2' }.map { id, cls, a -> [ id, a ] },
               remainder: true )
        .filter { id, a1, a2 -> a1 != null }
        .combine( Channel.fromList(arms) )
        .map { id, a1, a2, arm -> [ [id: id, cls: 'mhc1', arm: arm], [ mhc1: a1, mhc2: a2 ?: '' ] ] }

    // ---- the two arms -------------------------------------------------------------------------
    //
    // `remainder: true` on the rerank join so a table with no window FASTA still runs: it simply
    // has no wild type, which `wt_absent` carries, and which is the honest state rather than an
    // imputed one.
    ch_rerank = !('rerank' in arms) ? Channel.empty() :
        ch_epi.map { id, cls, f -> [ [id: id, cls: cls], f ] }
              .join( ch_win.map { id, cls, f -> [ [id: id, cls: cls], f ] }, remainder: true )
              .filter { meta, tsv, fa -> tsv != null }
              .map { meta, tsv, fa ->
                  [ meta + [arm: 'rerank'], tsv, fa ?: file("${moduleDir}/NO_FILE"), meta.cls ] }

    ch_denovo = !('denovo' in arms) ? Channel.empty() :
        ch_win.combine( ch_lists, by: [0, 1] )
              .map { id, cls, f, a -> [ [id: id, cls: cls, arm: 'denovo'], f, a, cls ] }
              // No allele list, nothing to score against. `predict` REQUIRES --alleles, and a
              // silently empty panel is the failure `mhcmatch alleles` exists to prevent -- so a
              // sample without one is dropped loudly here rather than scored against nothing.
              .filter { meta, f, a, cls ->
                  if( !a ) log.warn "no alleles for ${meta.id} ${cls}: skipping the de novo arm"
                  a as Boolean }

    MHCMATCH( ch_rerank, ch_denovo, ch_alleles )
}
