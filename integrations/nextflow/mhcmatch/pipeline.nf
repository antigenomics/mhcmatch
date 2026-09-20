#!/usr/bin/env nextflow
//
// mhcmatch, runnable from a samplesheet.
//
//   nextflow run pipeline.nf --input samplesheet.csv --outdir results --mode both \
//       --mhcmatch_vector_n0 8
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
//               `.genotype.json`, or one allele per line). Omitted for every row, and
//               `--alleles` / `--alleles_mhc2` must give a literal list instead.
//
// **A relative cell resolves against the SAMPLESHEET's own directory, not the launch directory.**
// A sheet is written beside the files it names and then run from wherever the work happens; the
// launch directory is not a property of the data and resolving against it makes the same sheet
// mean different things from two terminals.
//
// The assumed upstream is the community-standard stack: **nf-core/sarek -> VEP -> pVACtools** for
// the variants and the candidate table, **OptiType** (class I) and **arcasHLA / HLA-LA** (class II)
// for the typing. `pvacseq generate_protein_fasta` is what produces the peptide-window FASTA.

nextflow.enable.dsl = 2

include { MHCMATCH_ALLELES        } from './main.nf'
include { MHCMATCH_RERANK_ARM     } from './subworkflows/rerank.nf'
include { MHCMATCH_DENOVO_ARM     } from './subworkflows/denovo.nf'

// **Everything below the includes is a declaration, not a statement.** Nextflow 26.x strict syntax
// rejects a bare statement at script level -- "Statements cannot be mixed with script declarations"
// -- so `params.input = null` and `def sheetPath = { ... }` both fail to compile. The params this
// script owns are declared in ./nextflow.config beside every `params.mhcmatch_*`, and this helper is
// a FUNCTION rather than a closure assigned to a name.

//: One samplesheet cell -> a staged path, or null when the cell is empty.
//:
//: A relative cell is resolved against `base`, the samplesheet's own parent, for the reason in the
//: header. `checkIfExists` is on deliberately: Nextflow stages a `path` input by symlink and does
//: not check the target, so a mistyped cell otherwise becomes a DANGLING link that a `-stub-run`
//: never reads and reports success on, and the real run then dies inside a task on a bare basename.
//: A sheet is the one place where that typo is cheap to catch.
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

    def wantRerank = params.mode in ['rerank', 'both']
    def wantDenovo = params.mode in ['denovo', 'both']

    // ---- the samplesheet, validated per row ---------------------------------------------------
    ch_rows = Channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row ->
            // **`row['class']` rather than `row.class`**, for the one column whose name is also a
            // Groovy property. Measured on Groovy 5.0.8: a Map's property access *does* return the
            // entry, so `row.class` yields `mhc1` and both spellings work today -- this is not a
            // bug being worked around. The subscript form is used because it cannot be read, by a
            // reviewer or by a stricter parser, as the `getClass()` it looks like, and it is what
            // nf-core samplesheets use for exactly this collision.
            def id  = row.sample?.toString()?.trim()
            def cls = row['class']?.toString()?.trim()
            if( !id )
                error "samplesheet ${sheet.name}: a row has no `sample`. The columns are sample,class,candidates,windows,hla"
            if( !(cls in ['mhc1', 'mhc2']) )
                error "samplesheet ${sheet.name}, sample '${id}': `class` must be mhc1 or mhc2 (got '${cls ?: ''}')"
            def cand = sheetPath(base, row.candidates)
            def win  = sheetPath(base, row.windows)
            def hla  = sheetPath(base, row.hla)
            // A row naming neither input has nothing for either arm to read. Dropping it quietly is
            // how a cohort silently loses a donor; the sheet is the one artifact that can say so by
            // name, before any work starts.
            if( !cand && !win )
                error "samplesheet ${sheet.name}, sample '${id}' ${cls}: neither `candidates` nor " +
                      "`windows` is given, so there is nothing to score. `candidates` feeds the " +
                      "rerank arm, `windows` feeds the de novo arm and the rerank arm's --context; " +
                      "one of the two is required"
            [ id, cls, cand, win, hla ]
        }

    // **An empty input is the failure this whole module keeps meeting**, and it is not allowed to be
    // silent: a header-only sheet otherwise produces a run that does nothing and exits 0.
    ch_rows.count().subscribe { n ->
        if( n == 0 )
            error "samplesheet ${params.input} has a header and no data rows"
        else
            log.info "mhcmatch: ${n} (sample, class) row(s) from ${sheet.name}, mode=${params.mode}"
    }

    // One row per (sample, class), because every process writes `${meta.id}.<something>` and a
    // duplicate row runs the whole chain twice under ONE output name -- two tasks publishing the
    // same filename, which is the collision `nextflow.config`'s arm prefixes exist to prevent.
    //
    // One param and an index, NOT `{ key, hits -> }`: `map` destructures a tuple into a
    // multi-parameter closure, `subscribe` hands the item over WHOLE, and the two-param form
    // silently binds the entire pair to the first name.
    ch_rows.map { id, cls, cand, win, hla -> [ [id, cls], id ] }
           .groupTuple()
           .subscribe { hit ->
               if( (hit[1] as List).size() > 1 )
                   error "samplesheet ${sheet.name} names ${hit[0][0]} ${hit[0][1]} " +
                         "${(hit[1] as List).size()} times; one row per (sample, class)"
           }

    ch_epi  = ch_rows.filter { id, cls, cand, win, hla -> cand != null }
                     .map    { id, cls, cand, win, hla -> [ id, cls, cand ] }
    ch_win  = ch_rows.filter { id, cls, cand, win, hla -> win != null }
                     .map    { id, cls, cand, win, hla -> [ id, cls, win ] }
    ch_keys = ch_rows.map { id, cls, cand, win, hla -> [ id, cls ] }.unique()

    // A mode with no usable row is a typo in the sheet or in `--mode`, not an empty result.
    if( wantRerank )
        ch_epi.count().subscribe { n ->
            if( n == 0 )
                error "--mode ${params.mode} needs a `candidates` table and no row of ${sheet.name} names one" }
    if( wantDenovo )
        ch_win.count().subscribe { n ->
            if( n == 0 )
                error "--mode ${params.mode} needs a `windows` FASTA and no row of ${sheet.name} names one" }

    // A typing file is a property of the DONOR, not of the class, so both of a sample's rows name
    // the same one and it is resolved once per sample. Two different files under one sample is a
    // sheet-authoring mistake that would otherwise run `mhcmatch alleles` twice for one (sample,
    // class) key and publish both under one name.
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
                           [ id, distinct[0] ]
                       }

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
    else {
        // No typing file anywhere and no literal is a run with no panel at all, which `predict`
        // refuses by name and every other path answers with an empty panel and exit 0.
        ch_typing.count().subscribe { n ->
            if( n == 0 )
                error "no row of ${sheet.name} names an `hla` file, so there is nothing to resolve " +
                      "a panel from. Fill the `hla` column, or pass --alleles / --alleles_mhc2 with " +
                      "a literal list for every sample (the mouse case: an inbred line's haplotype " +
                      "is a property of the line)" }

        MHCMATCH_ALLELES( ch_keys.combine( ch_typing, by: 0 )
                                 .map { id, cls, f -> [ [id: id, cls: cls], f, cls ] } )
        // **The rows WITHOUT a typing file are mixed back in, carrying an empty list.** The join
        // above is an inner one, so on its own it drops those keys entirely -- and a key missing
        // from this channel vanishes from the de novo arm without a word, which is the silent
        // version of the loud warning at that arm's filter. The samplesheet says per row which
        // case a sample is in, so neither has to be inferred from an absence.
        // **`groupTuple` then prefer the non-empty, because a key can arrive twice.** `ch_typing`
        // is keyed on `id` alone, so a sample whose mhc1 row names an `hla` file and whose mhc2 row
        // leaves it blank gets a resolved list from the combine AND an empty one from the mix.
        // Before the group, the rerank arm's `join` paired whichever buffered first -- and the
        // empty one is emitted immediately while the resolved one waits on a process -- so
        // `self_help` was computed against an empty class-II panel with nothing said. Grouping
        // makes the result independent of arrival order, which is the property that was missing.
        ch_alleles = MHCMATCH_ALLELES.out.alleles
                        .map { meta, cls, f -> [ meta.id, cls, f.text.trim() ] }
                        .mix( ch_rows.filter { id, cls, cand, win, hla -> hla == null }
                                     .map    { id, cls, cand, win, hla -> [ id, cls, '' ] } )
                        .groupTuple( by: [0, 1] )
                        .map { id, cls, lists -> [ id, cls, lists.find { it } ?: '' ] }
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
            // The donor's class-I list AND their class-II list, as one value. MHCMATCH_CASSETTE
            // takes `[mhc1:, mhc2:]` and computes `self_help` from the second -- whether a unit's
            // CD8 epitope has CD4 help from the SAME unit, which is what the cassette map is for.
            // A per-donor class-II list has no other way in: a sixth element on that process's
            // input tuple would break every pipeline that `include`s it.
            //
            // `remainder: true` keeps a donor with no class-II input, who simply gets no
            // `self_help`; the filter drops the mirror case, an id with class II and no class I,
            // which has nothing to build a cassette from.
            ch_alleles.filter { id, cls, a -> cls == 'mhc1' }.map { id, cls, a -> [ id, a ] }
                      .join( ch_alleles.filter { id, cls, a -> cls == 'mhc2' }
                                       .map { id, cls, a -> [ id, a ] }, remainder: true )
                      .filter { id, a1, a2 -> a1 != null }
                      .map { id, a1, a2 -> [ [id: id, cls: 'mhc1'], [ mhc1: a1, mhc2: a2 ?: '' ] ] }
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
