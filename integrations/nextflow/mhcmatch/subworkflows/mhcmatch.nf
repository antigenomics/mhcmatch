// The two arms and the tail they share. ONE workflow, because the tail is the same five processes
// over the same donor either way.
//
//   rerank:  your table  ─► RERANK ─► *.epitopes.mhcmatch.tsv ─┐
//   denovo:  windows.fa  ─► PREDICT ─► scored.csv + native.tsv │
//                        └► RANK ───► *.mhcmatch.ranked.tsv ───┤
//                                                              │  (class I only, below here)
//              ┌───────────────────────────────────────────────┘
//              ├─► NEOAG           ─► neoag.tsv
//              ├─► MIMICRY         ─► mimicry.tsv        (off by default)
//              └─► CASSETTE_SELECT ─► units.tsv
//                    └─► CASSETTE  ─► .faa / .fna / map
//                          └─► CASSETTE_SCORE            ONE per arm, over EVERY donor
//
// **The arm rides in `meta.arm`, not in an alias.** A DSL2 process may be invoked once per run, so
// two arms used to mean a second `_DN`-suffixed include of every process, a duplicate copy of this
// file, and six `withName:` selectors in nextflow.config whose job was to tell the copies apart --
// a selector spelled as the bare name sizes one arm and silently misses the other. One instance
// fed by two channels has neither problem: the arm is a field, `armPrefix` reads it for the output
// name, and `CASSETTE_SCORE` groups on it so each arm still gets its own cohort calibration.
//
// **Everything below the fork is CLASS I ONLY, by design and not by omission.** Prior evidence and
// safety are built on a CD8 mechanism -- a minimal epitope close enough to a confirmed neoantigen
// that one clonotype could see both, and a register that IS an essential-tissue self peptide,
// killing the cell presenting it. Neither becomes a class-II question by widening the length
// range: CD4 self-reactivity runs through help, hypersensitivity and allergy, which has different
// thresholds and none of them measured here. `select` also spends per-allotype capacity, and the
// class-II locus call is not good enough to spend it on (mhcmatch and ISP agree on the presenting
// locus for 52.7% of class-II rows against 78.1% for class I). See docs/safety.rst.

include { MHCMATCH_PREDICT         } from '../main.nf'
include { MHCMATCH_RANK            } from '../main.nf'
include { MHCMATCH_RERANK          } from '../main.nf'
include { MHCMATCH_NEOAG           } from '../main.nf'
include { MHCMATCH_MIMICRY         } from '../main.nf'
include { MHCMATCH_CASSETTE_SELECT } from '../main.nf'
include { MHCMATCH_CASSETTE        } from '../main.nf'
include { MHCMATCH_CASSETTE_SCORE  } from '../main.nf'

//: The class-I list, whichever of the two shapes the allele value has -- a bare String, or
//: `[mhc1: '...', mhc2: '...']` when the caller also has the donor's class-II allotypes. Only
//: MHCMATCH_CASSETTE reads the class-II half, to compute `self_help` for the cassette map.
def mhc1Of(a) {
    a instanceof Map ? (a.mhc1 ?: '') : (a ?: '')
}

workflow MHCMATCH {

    take:
    ch_rerank    // [ val(meta), path(table),  path(context|NO_FILE), val(cls) ]  -- may be empty
    ch_denovo    // [ val(meta), path(fasta),  val(alleles),          val(cls) ]  -- may be empty
    ch_alleles   // [ val(meta), val(alleles) ] for the class-I meta of each donor. A String, or
                 // `[mhc1: '...', mhc2: '...']` to give the cassette map the class-II list too.

    main:
    // **Forgetting the `includeConfig` line is the failure this catches.** Nextflow auto-loads the
    // config beside the ENTRY script, so `nextflow run pipeline.nf` gets it for free and an outside
    // pipeline that only `include`s this file does not. Every `params.mhcmatch_*` is then undefined,
    // which Nextflow reports as one WARN among many and evaluates as null -- and `isOn(null)` is
    // FALSE, so the cassette safety screen reads as *off* and MHCMATCH_CASSETTE builds a construct
    // with no screen at all. A WARN is not enough warning for that. `containsKey`, not truth: a
    // caller who deliberately set `false` has made a choice and is left alone.
    if( !params.containsKey('mhcmatch_cassette_screen') )
        error "mhcmatch: params.mhcmatch_cassette_screen is not defined, so the module's config " +
              "was never loaded -- and the cassette safety screen would run OFF without saying " +
              "so. Add `includeConfig '<path>/integrations/nextflow/mhcmatch/nextflow.config'` to " +
              "your pipeline's config (after your own params), or define the mhcmatch_* params."

    MHCMATCH_RERANK( ch_rerank )
    MHCMATCH_PREDICT( ch_denovo )
    MHCMATCH_RANK( ch_denovo )

    // The two arms' scored tables, one shape, tagged by `meta.arm`. Everything below is arm-blind.
    ch_pool = MHCMATCH_RERANK.out.reranked
        .mix( MHCMATCH_RANK.out.ranked )
        .filter { meta, cls, tsv -> cls == 'mhc1' }
        .map    { meta, cls, tsv -> [ meta, tsv ] }

    MHCMATCH_NEOAG( ch_pool.map { meta, tsv -> [ meta, tsv, 'mhc1' ] } )
    // **Mimicry is an ANNOTATION step and ships off.** It says what a candidate resembles; it does
    // not feed the ranking, because `rank`'s corpus channels are a `corpus_spectrum` contraction
    // rather than a neighbour search and build no index at all. Turning it on costs a whole-
    // proteome reference index -- 65.0 s warm-cached, ~194 s per task cold, and in a fan-out every
    // task misses the cache at once and builds it simultaneously. That was the dominant stage of
    // the whole pipeline.
    MHCMATCH_MIMICRY( params.mhcmatch_mimicry?.toString()?.toLowerCase() in
                          [null, 'false', '0', 'no', 'null', '']
                      ? Channel.empty()
                      : ch_pool.map { meta, tsv -> [ meta, tsv, 'mhc1' ] } )

    // `remainder: true` so a sample with no allele list still reaches the selector -- it loses the
    // allotype channel and its coverage denominator, and says so. The filter drops the mirror case
    // the same flag emits: an allele entry that matched no pool, which is a real shape when a
    // donor has a window FASTA and no candidate table.
    ch_scored = ch_pool.join( ch_alleles, remainder: true )
                       .filter { meta, tsv, alleles -> tsv != null }

    MHCMATCH_CASSETTE_SELECT( ch_scored.map { meta, tsv, a -> [ meta, tsv, mhc1Of(a) ] } )

    // The window FASTA, carried to `--context`: `rank`/`rerank` emit MINIMAL epitopes and a unit is
    // the long (~27 aa) window around the mutation, so neither side alone can build one. A row that
    // genuinely names no window still arrives as the NO_FILE sentinel, so `--unit-column` remains
    // the other source and MHCMATCH_CASSETTE's refusal survives for a sample with neither.
    //
    // `alleles` goes in WHOLE here where the selector got only the class-I half: this is the one
    // process that reads the class-II list.
    ch_ctx = ch_rerank.filter { meta, tsv, ctx, cls -> cls == 'mhc1' }
                      .map    { meta, tsv, ctx, cls -> [ meta, ctx ] }
             .mix( ch_denovo.filter { meta, fa, a, cls -> cls == 'mhc1' }
                            .map    { meta, fa, a, cls -> [ meta, fa ] } )

    MHCMATCH_CASSETTE(
        MHCMATCH_CASSETTE_SELECT.out.units
            .join( ch_alleles, remainder: true )
            .filter { meta, units, alleles -> units != null }
            .join( ch_ctx, remainder: true )
            // `moduleDir` and not `projectDir`: projectDir is the ENTRY script's directory, so an
            // integrator including this subworkflow resolves the sentinel against THEIR repo root,
            // where it does not exist.
            .map { meta, units, alleles, ctx ->
                [ meta, units, ctx ?: file("${moduleDir}/../NO_FILE"), alleles ?: '', 'mhc1' ] }
    )

    // ONE calibration per arm over every donor in it, which is why this collects -- see the note on
    // the process. It takes the **units** table and not the `.cassette.tsv` report: `cassette
    // score` wants one row per manufactured unit with a peptide and a score, and the report is
    // long-form with neither. Joined on CASSETTE so the score still waits for assembly -- a
    // cassette that failed its safety screen should not be scored as if it shipped.
    MHCMATCH_CASSETTE_SCORE(
        MHCMATCH_CASSETTE_SELECT.out.units
            .join( MHCMATCH_CASSETTE.out.report )
            .map { meta, units, report -> [ meta.arm, units ] }
            .groupTuple()
            .join( ch_pool.map { meta, tsv -> [ meta.arm, tsv ] }.groupTuple() )
    )

    emit:
    scored     = MHCMATCH_PREDICT.out.scored          // [ meta, cls, *.mhcmatch.scored.csv ]
    native_tsv = MHCMATCH_PREDICT.out.native_tsv      // [ meta, cls, *.mhcmatch.native.tsv ]
    ranked     = MHCMATCH_RANK.out.ranked             // [ meta, cls, *.mhcmatch.ranked.tsv ]
    reranked   = MHCMATCH_RERANK.out.reranked         // [ meta, cls, *.epitopes.mhcmatch.tsv ]
    neoag      = MHCMATCH_NEOAG.out.neoag
    mimicry    = MHCMATCH_MIMICRY.out.mimicry
    units      = MHCMATCH_CASSETTE_SELECT.out.units   // [ meta, *.vaccine.units.tsv ]
    cassette   = MHCMATCH_CASSETTE.out.protein        // [ meta, *.cassette.faa ]
    cds        = MHCMATCH_CASSETTE.out.cds            // [ meta, *.cassette.fna ]
    report     = MHCMATCH_CASSETTE.out.report
    map        = MHCMATCH_CASSETTE.out.map
    score      = MHCMATCH_CASSETTE_SCORE.out.score    // [ arm, cohort.<arm>.cassette_score.tsv ]
}
