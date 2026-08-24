// End-to-end mhcmatch: variant windows in, ranked candidates and a screened cassette out.
//
// Chains the five processes in ./main.nf. Take it as written or as a template -- the useful part is
// which output feeds which input, because that is the thing a README sentence gets wrong.
//
//   windows.fasta ─► PREDICT ─► native.tsv                        both classes
//                 └► RANK ────► ranked.tsv ─┬─► NEOAG   ─► neoag.tsv     class I only
//                                           ├─► MIMICRY ─► mimicry.tsv   class I only
//                                           └─► CASSETTE─► cassette      class I only
//                                                        └─► CASSETTE_SCORE (all donors at once)
//
// PREDICT and RANK serve both classes. NEOAG, MIMICRY and CASSETTE are CD8-only by design, not
// by omission -- see the comment at the filter and docs/safety.rst.
//
// CASSETTE_SCORE is the one step that is NOT per sample: it waits for every donor so the
// calibration offset is fitted once over the run. See its comment in ../main.nf.
//
// CASSETTE takes `ranked.tsv` as its candidate table AND the original `windows.fasta` as `--context`:
// `rank` emits minimal epitopes and a unit is the long window around the mutation, so neither side
// alone can build one (see mhcmatch.vector.units_from_context). Injecting the minimal epitope is not
// a smaller version of the right thing, it is the tolerising configuration.

include { MHCMATCH_PREDICT  } from '../main.nf'
include { MHCMATCH_RANK     } from '../main.nf'
include { MHCMATCH_NEOAG    } from '../main.nf'
include { MHCMATCH_MIMICRY  } from '../main.nf'
include { MHCMATCH_CASSETTE       } from '../main.nf'
include { MHCMATCH_CASSETTE_SCORE } from '../main.nf'

workflow MHCMATCH {

    take:
    ch_windows        // [ val(meta), path(fasta), val(alleles), val(cls) ]

    main:
    ch_versions = Channel.empty()

    MHCMATCH_PREDICT( ch_windows )
    ch_versions = ch_versions.mix( MHCMATCH_PREDICT.out.versions.first() )

    MHCMATCH_RANK( ch_windows )
    ch_versions = ch_versions.mix( MHCMATCH_RANK.out.versions.first() )

    // `neoag` and `mimicry` read a table with a `peptide` column and carry every other column
    // through, so the ranked table goes in unchanged and comes back annotated.
    ch_ranked = MHCMATCH_RANK.out.ranked                       // [ meta, cls, ranked.tsv ]

    // CLASS I ONLY, for all three of NEOAG, MIMICRY and CASSETTE. Prior evidence and safety are built
    // on a CD8 mechanism -- a minimal epitope close enough to a confirmed neoantigen that one
    // clonotype could see both, and a register that IS an essential-tissue self peptide killing the
    // cell presenting it. Neither becomes a class-II question by widening the length range: CD4
    // self-reactivity runs through help, hypersensitivity and allergy, which has different
    // thresholds and none of them measured here. Running these on class II would produce columns
    // that look like answers and are not. See docs/safety.rst.
    ch_pep    = ch_ranked
        .filter { meta, cls, tsv -> cls == 'mhc1' }
        .map    { meta, cls, tsv -> [ meta, tsv, cls ] }

    MHCMATCH_NEOAG( ch_pep )
    MHCMATCH_MIMICRY( ch_pep )
    ch_versions = ch_versions.mix( MHCMATCH_NEOAG.out.versions.first() )
    ch_versions = ch_versions.mix( MHCMATCH_MIMICRY.out.versions.first() )

    // Cassette assembly is class-I only for the reason above, and for a second one: `select`
    // spends per-allotype capacity, and the class-II locus call is not good enough to spend it on
    // (mhcmatch and ISP agree on the presenting locus for 52.7% of class-II rows against 78.1% for
    // class I).
    ch_vector = ch_ranked
        .filter { meta, cls, tsv -> cls == 'mhc1' }
        .join( ch_windows.map { meta, fa, alleles, cls -> [ meta, fa, alleles ] } )
        .map { meta, cls, ranked, fa, alleles -> [ meta, ranked, fa, alleles, cls ] }

    MHCMATCH_CASSETTE( ch_vector )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE.out.versions.first() )

    // ONE calibration for the whole run, which is why this is `.collect()` and not a per-sample
    // call. `rank` anchors `p_response` on the batch it is handed, so scoring each donor alone
    // makes every donor's mean the declared prevalence and no two donors comparable. Collecting
    // first is the fix, and it is the only place in this subworkflow where a process deliberately
    // waits for every sample.
    MHCMATCH_CASSETTE_SCORE(
        MHCMATCH_CASSETTE.out.report.map { meta, tsv -> tsv }.collect(),
        ch_ranked.filter { meta, cls, tsv -> cls == 'mhc1' }.map { meta, cls, tsv -> tsv }.collect()
    )
    ch_versions = ch_versions.mix( MHCMATCH_CASSETTE_SCORE.out.versions )

    emit:
    scored   = MHCMATCH_PREDICT.out.scored     // [ meta, cls, *.mhcmatch.scored.csv ]
    native_tsv = MHCMATCH_PREDICT.out.native_tsv   // [ meta, cls, *.mhcmatch.native.tsv ]
    ranked   = MHCMATCH_RANK.out.ranked        // [ meta, cls, *.mhcmatch.ranked.tsv ]
    neoag    = MHCMATCH_NEOAG.out.neoag        // [ meta, cls, *.mhcmatch.neoag.tsv ]
    mimicry  = MHCMATCH_MIMICRY.out.mimicry    // [ meta, cls, *.mhcmatch.mimicry.tsv ]
    cassette = MHCMATCH_CASSETTE.out.protein   // [ meta, *.cassette.faa ]
    cds      = MHCMATCH_CASSETTE.out.cds       // [ meta, *.cassette.fna ]
    report   = MHCMATCH_CASSETTE.out.report    // [ meta, *.cassette.tsv ]
    score    = MHCMATCH_CASSETTE_SCORE.out.score  // cohort.cassette_score.tsv, ONE per run
    versions = ch_versions
}
