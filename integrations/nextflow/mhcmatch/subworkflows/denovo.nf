// The DE NOVO arm: a mutation-window FASTA in, an epitope table built entirely by mhcmatch out,
// then the cassette.
//
//   windows.fasta ─► PREDICT ─► scored.csv + native.tsv
//                 └► RANK ────► ranked.tsv ─┬─► NEOAG   ─► neoag.tsv        class I
//                                           ├─► MIMICRY ─► mimicry.tsv      class I
//                                           └─► CASSETTE_SELECT ─► units.tsv
//                                                 └─► CASSETTE ─► .faa / .fna / map
//                                                       └─► CASSETTE_SCORE (all donors at once)
//
// Where ./mhcmatch.nf differs: that subworkflow hands `ranked.tsv` straight to CASSETTE, whose
// `--n0` stopping rule decides how many units the recipient's allotypes can carry. This one puts
// CASSETTE_SELECT in front, so the size is `-k` -- what a trial that has committed to a construct
// size actually needs. Both are kept because both are real questions; ./mhcmatch.nf is unchanged so
// an existing `include` of it still works.
//
// CASSETTE takes `ranked.tsv` AND the original `windows.fasta` as `--context`: `rank` emits minimal
// epitopes and a unit is the long window around the mutation, so neither side alone can build one.

// **Aliased, and that is not cosmetic.** A DSL2 process may be *invoked* once per run; including
// it in two subworkflows and running both raises "Process 'X' has been already used". The rerank
// arm takes the plain names, so an existing config selector keeps matching it; this arm takes a
// `_DN` suffix, which also makes the two arms tellable apart in the trace under `--mode both`. The
// `withName:` selectors in ../nextflow.config are written to match either.
include { MHCMATCH_PREDICT                                     } from '../main.nf'
include { MHCMATCH_RANK                                        } from '../main.nf'
include { MHCMATCH_NEOAG           as MHCMATCH_NEOAG_DN        } from '../main.nf'
include { MHCMATCH_MIMICRY         as MHCMATCH_MIMICRY_DN      } from '../main.nf'
include { MHCMATCH_CASSETTE_SELECT as MHCMATCH_CASSETTE_SELECT_DN } from '../main.nf'
include { MHCMATCH_CASSETTE        as MHCMATCH_CASSETTE_DN     } from '../main.nf'
include { MHCMATCH_CASSETTE_SCORE  as MHCMATCH_CASSETTE_SCORE_DN } from '../main.nf'

workflow MHCMATCH_DENOVO_ARM {

    take:
    ch_windows        // [ val(meta), path(fasta), val(alleles), val(cls) ]

    main:
    ch_versions = Channel.empty()

    MHCMATCH_PREDICT( ch_windows )
    MHCMATCH_RANK( ch_windows )
    ch_versions = ch_versions.mix( MHCMATCH_PREDICT.out.versions.first() )
    ch_versions = ch_versions.mix( MHCMATCH_RANK.out.versions.first() )

    // CLASS I ONLY below -- see the comment at the same filter in ./mhcmatch.nf and docs/safety.rst.
    ch_mhc1 = MHCMATCH_RANK.out.ranked.filter { meta, cls, tsv -> cls == 'mhc1' }

    MHCMATCH_NEOAG_DN( ch_mhc1.map { meta, cls, tsv -> [ meta, tsv, cls ] } )
    MHCMATCH_MIMICRY_DN( ch_mhc1.map { meta, cls, tsv -> [ meta, tsv, cls ] } )
    ch_versions = ch_versions.mix( MHCMATCH_NEOAG_DN.out.versions.first() )
    ch_versions = ch_versions.mix( MHCMATCH_MIMICRY_DN.out.versions.first() )

    // The class-I window FASTA and allele list, keyed by sample, for the two cassette steps -- with
    // the SAME donor's class-II list folded in beside it, because this arm already has it: a de novo
    // run over both classes carries a `cls == 'mhc2'` row for the same sample, holding exactly the
    // allotypes the cassette map needs for `self_help`. MHCMATCH_CASSETTE takes the pair as a Map
    // (../main.nf), which is what keeps a per-donor list off the input tuple.
    //
    // Joined on `meta.id`, not on `meta`, because the two rows' metas differ by `cls`. And
    // `remainder: true`: a donor with no class-II windows still gets a cassette, with the map
    // falling back to `params.mhcmatch_vector_map_alleles_mhc2` and saying so if that is unset too.
    ch_a2  = ch_windows.filter { meta, fa, alleles, cls -> cls == 'mhc2' }
                       .map    { meta, fa, alleles, cls -> [ meta.id, alleles ] }
    ch_ctx = ch_windows.filter { meta, fa, alleles, cls -> cls == 'mhc1' }
                       .map    { meta, fa, alleles, cls -> [ meta.id, meta, fa, alleles ] }
                       .join( ch_a2, remainder: true )
                       .filter { id, meta, fa, a1, a2 -> meta != null }
                       .map    { id, meta, fa, a1, a2 ->
                           [ meta, fa, [ mhc1: a1, mhc2: a2 ?: '' ] ] }

    // The selector prices class-I allotype coverage, so it takes the class-I half only.
    MHCMATCH_CASSETTE_SELECT_DN(
        ch_mhc1.map { meta, cls, tsv -> [ meta, tsv ] }
               .join( ch_ctx.map { meta, fa, alleles -> [ meta, alleles.mhc1 ] } )
    )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE_SELECT_DN.out.versions.first() )

    MHCMATCH_CASSETTE_DN(
        MHCMATCH_CASSETTE_SELECT_DN.out.units
            .join( ch_ctx )
            .map { meta, units, fa, alleles -> [ meta, units, fa, alleles, 'mhc1' ] }
    )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE_DN.out.versions.first() )

    // The **units** table, not the `.cassette.tsv` report -- see the same call in ./rerank.nf.
    MHCMATCH_CASSETTE_SCORE_DN(
        MHCMATCH_CASSETTE_SELECT_DN.out.units
            .join( MHCMATCH_CASSETTE_DN.out.report )
            .map { meta, units, report -> units }.collect(),
        ch_mhc1.map { meta, cls, tsv -> tsv }.collect()
    )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE_SCORE_DN.out.versions )

    emit:
    scored     = MHCMATCH_PREDICT.out.scored
    native_tsv = MHCMATCH_PREDICT.out.native_tsv
    ranked     = MHCMATCH_RANK.out.ranked
    neoag      = MHCMATCH_NEOAG_DN.out.neoag
    mimicry    = MHCMATCH_MIMICRY_DN.out.mimicry
    units      = MHCMATCH_CASSETTE_SELECT_DN.out.units
    cassette   = MHCMATCH_CASSETTE_DN.out.protein
    cds        = MHCMATCH_CASSETTE_DN.out.cds
    report     = MHCMATCH_CASSETTE_DN.out.report
    score      = MHCMATCH_CASSETTE_SCORE_DN.out.score
    versions   = ch_versions
}
