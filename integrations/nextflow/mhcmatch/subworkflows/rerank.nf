// The RERANK arm: a caller's own candidate table in, the same table out with this model appended
// and re-ordered by it, then the cassette.
//
//   epitopes.tsv ─► RERANK ─► *.epitopes.mhcmatch.tsv ─┬─► NEOAG    ─► neoag.tsv      class I
//   (+ windows.fasta as --context)                     ├─► MIMICRY  ─► mimicry.tsv    class I
//                                                      └─► CASSETTE_SELECT ─► units.tsv
//                                                            └─► CASSETTE ─► .faa / .fna / map
//                                                                  └─► CASSETTE_SCORE (all donors)
//
// The whole point of this arm is that **the caller's columns survive**. `rank --passthrough` emits
// them verbatim and in the caller's order, ahead of ours, so the deliverable is their table with a
// column block added -- not a different table they then have to join back.
//
// CASSETTE reads `epitope_context` as its unit (`--unit-column`), because the rerank arm has no
// window FASTA to rebuild a unit from and the table already carries one at 27 aa. That is the same
// object `--context` produces on the de novo arm, arrived at from the other side: a minimal epitope
// loads onto any cell without costimulation and is the tolerising configuration, so neither path
// is allowed to inject one.

include { MHCMATCH_RERANK          } from '../main.nf'
include { MHCMATCH_NEOAG           } from '../main.nf'
include { MHCMATCH_MIMICRY         } from '../main.nf'
include { MHCMATCH_CASSETTE_SELECT } from '../main.nf'
include { MHCMATCH_CASSETTE        } from '../main.nf'
include { MHCMATCH_CASSETTE_SCORE  } from '../main.nf'

workflow MHCMATCH_RERANK_ARM {

    take:
    ch_tables         // [ val(meta), path(table), path(context|NO_FILE), val(cls) ]
    ch_alleles        // [ val(meta), val(alleles) ] -- the donor's DISTINCT class-I allotypes

    main:
    ch_versions = Channel.empty()

    MHCMATCH_RERANK( ch_tables )
    ch_versions = ch_versions.mix( MHCMATCH_RERANK.out.versions.first() )

    // CLASS I ONLY for NEOAG, MIMICRY and the cassette -- by design, not omission. Prior evidence
    // and safety are built on a CD8 mechanism; CD4 self-reactivity runs through help,
    // hypersensitivity and allergy, which has different thresholds and none of them measured here.
    // See docs/safety.rst and the same filter in ./mhcmatch.nf.
    ch_mhc1 = MHCMATCH_RERANK.out.reranked.filter { meta, cls, tsv -> cls == 'mhc1' }
    ch_pep  = ch_mhc1.map { meta, cls, tsv -> [ meta, tsv, cls ] }

    MHCMATCH_NEOAG( ch_pep )
    MHCMATCH_MIMICRY( ch_pep )
    ch_versions = ch_versions.mix( MHCMATCH_NEOAG.out.versions.first() )
    ch_versions = ch_versions.mix( MHCMATCH_MIMICRY.out.versions.first() )

    ch_pool = ch_mhc1.map { meta, cls, tsv -> [ meta, tsv ] }

    // `remainder: true` so a sample with no allele list still reaches the selector (it loses the
    // allotype channel and its coverage denominator, and says so) -- but the SAME flag also emits
    // an allele entry that matched no pool, which is a real shape under `--mode both` over a mixed
    // directory: a donor with a window FASTA and no candidate table has a class-I allele list and
    // nothing to rerank. Handing that through gives the process a null path, so it is filtered.
    MHCMATCH_CASSETTE_SELECT( ch_pool.join( ch_alleles, remainder: true )
                                     .filter { meta, tsv, alleles -> tsv != null }
                                     .map { meta, tsv, alleles -> [ meta, tsv, alleles ?: '' ] } )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE_SELECT.out.versions.first() )

    // The chosen units carry the caller's columns through `select --passthrough`, so the long
    // window is still there for CASSETTE to build from. `NO_FILE` for context: there is none.
    MHCMATCH_CASSETTE(
        MHCMATCH_CASSETTE_SELECT.out.units
            .join( ch_alleles, remainder: true )
            .filter { meta, units, alleles -> units != null }
            .map { meta, units, alleles ->
                [ meta, units, file("${projectDir}/NO_FILE"), alleles ?: '', 'mhc1' ] }
    )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE.out.versions.first() )

    // ONE calibration for the whole run, which is why this collects. `rank` anchors `p_response`
    // on the batch it is handed, so a per-donor call makes every donor's mean the declared
    // prevalence and no two donors comparable -- see the comment on the process in ../main.nf.
    //
    // It takes the **units** table and not the `.cassette.tsv` report: `cassette score` wants one
    // row per manufactured unit with a peptide and a score, and the report is long-form with
    // neither. Joined on CASSETTE so the score still waits for assembly -- a cassette that failed
    // its safety screen should not be scored as if it shipped.
    MHCMATCH_CASSETTE_SCORE(
        MHCMATCH_CASSETTE_SELECT.out.units
            .join( MHCMATCH_CASSETTE.out.report )
            .map { meta, units, report -> units }.collect(),
        ch_pool.map { meta, tsv -> tsv }.collect()
    )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE_SCORE.out.versions )

    emit:
    reranked = MHCMATCH_RERANK.out.reranked            // [ meta, cls, *.epitopes.mhcmatch.tsv ]
    neoag    = MHCMATCH_NEOAG.out.neoag
    mimicry  = MHCMATCH_MIMICRY.out.mimicry
    units    = MHCMATCH_CASSETTE_SELECT.out.units      // [ meta, *.vaccine.units.tsv ]
    cassette = MHCMATCH_CASSETTE.out.protein           // [ meta, *.cassette.faa ]
    cds      = MHCMATCH_CASSETTE.out.cds               // [ meta, *.cassette.fna ]
    report   = MHCMATCH_CASSETTE.out.report
    score    = MHCMATCH_CASSETTE_SCORE.out.score       // ONE per run
    versions = ch_versions
}
