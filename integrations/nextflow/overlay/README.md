# mhcmatch as an overlay on a neoantigen pipeline you already run

For a host pipeline that already does variant calling, HLA typing and expression quantification and
reaches a candidate table with those in hand — the **nf-core/sarek → VEP → pVACtools** shape, with
**OptiType** and **arcasHLA / HLA-LA** for the typing.

It attaches at named seams and **contributes zero new scoring processes**: everything is the module
in `../mhcmatch/`. The only two things this directory adds are `MHCMATCH_PREFLIGHT`, which computes
nothing and exists to refuse to start, and `MHCMATCH_RESCORE`, which writes the one column a host
selector reads.

## Why an overlay rather than copying the processes in

Because a copy is a fork, and a fork drifts. One hand-wired integration we reviewed had drifted five
ways within a few releases:

| drift | consequence |
|---|---|
| `-stub-run` hard-typed an 18-column header | stubs no longer matched the real 57-column table |
| it emitted `native_table` where ours emits `native_tsv` | a rename that only shows up when a downstream `include` is added |
| its image was pinned to mhcmatch **1.0.1** | sixteen releases behind, and nothing in a run says so |
| its `rank_threshold` defaulted to **2.0** | the weak class-I cut applied at the predictor, throwing away candidates the later aggregation could still have ranked. The same number is the **strong** cut for class II: measured on one window pair against `DRB1*15:01`, **0 of 56** scored pairs survive it |
| its class-II bridge called `mhcmatch predict` **once per allele**, in a Python loop | `predict` takes the whole allotype list and batches internally; N sequential calls also reload the model and rebuild the calibration N times |

An alias cannot drift. A sixth divergence will.

## Wiring it in

**Level 0 needs no edit to your repository at all** — one `-c`, plus your own predictor switch if
you have one, since a config file cannot reach inside your `if`:

```bash
nextflow run . -c /path/to/mhcmatch/integrations/nextflow/overlay/overlay.config \
    <your predictor switch> mhcmatch
```

That pins the image every mhcmatch-invoking process runs in and sets `--rank-threshold none`, so
filtering happens where you can see what was filtered. It changes nothing else, because
`--mhcmatch_overlay_mode` defaults to `off`.

For a pipeline with its own config, name **both** files, in this order, after your own params:

```groovy
includeConfig '<path>/integrations/nextflow/mhcmatch/nextflow.config'
includeConfig '<path>/integrations/nextflow/overlay/overlay.config'
```

Then include the workflow and wire your two seams:

```groovy
include { MHCMATCH_OVERLAY } from '<path>/integrations/nextflow/overlay/overlay.nf'

MHCMATCH_OVERLAY( ch_candidates, ch_windows, ch_alleles )

// the SELECTION seam: your selector reads this instead of your own table
ch_for_selector = MHCMATCH_OVERLAY.out.scored
// the CONSTRUCT seam: your codon optimiser reads ours only when we built one
ch_for_optimiser = params.mhcmatch_overlay_mode == 'cassette'
                   ? MHCMATCH_OVERLAY.out.vaccine : YOUR_CONSTRUCT.out.fasta
```

**A `withName` selector matches the alias.** If you `include { MHCMATCH_PREDICT as MHCMATCH_MHCI }`
— the ordinary way to serve two classes from one process — a selector naming `MHCMATCH_PREDICT`
binds to nothing and says nothing. So the image pin uses a deliberately broad literal,
`.*MHCMATCH_[A-Z0-9_]*$`; pinning the image on a process that never calls mhcmatch costs nothing.

## The five modes

| `--mhcmatch_overlay_mode` | `scored` emits | `vaccine` emits | what changes for you |
|---|---|---|---|
| `off` *(default)* | your table, verbatim | empty | nothing. Preflight only — safe to merge before anyone has decided |
| `annotate` | your table, verbatim | empty | a reranked table published **beside** yours. A no-op on every deliverable |
| `rerank` | your table with `score` rewritten and `score_host` carrying yours | empty | your selector, your linkers, your length budget, **our ordering** |
| `denovo` | your table, verbatim | empty | a third table from our de novo arm, orthogonal to yours |
| `cassette` | your table, verbatim | our `.cassette.faa` | our selector replaces yours. **Your own construct step still runs** — the ternary above just sends the optimiser a different FASTA |

`annotate` and `rerank` differ by exactly one process and that is the point of the mode: a host
selector sorts on a column whose name it fixes, so publishing `mm_score` beside it changes nothing
at all. `MHCMATCH_RESCORE` is what writes `score`.

## What `rerank` mode preserves

Your own scorer builds

```
score = base(affinity, agretopicity, expression) × driver × recurrence × hla_freq × hla_loss
```

and only `base` is a question about peptides. The multiplicative factors are **clinical** and this
model does not contain them: a driver boost of 100.0 and an HLA-loss weight of 0.0 for an allele the
tumour has lost are both ordinary settings, and replacing the product would put an epitope
restricted to a deleted allele back into a construct. So the base is replaced and the factors are
carried:

```
score := W × minmax(mm_score) × driver(driver_class) × hla_weight × recurrence_weight
         × hla_frequency_weight
```

`W` is your own total base weight, so the numbers stay in the range your operators read. Each factor
is used only if its column is present. Your value is kept as `score_host` — nothing is destroyed,
and the two are one `paste` apart.

Every host-specific name is a parameter with a default (`overlay.config`), so the process carries no
knowledge of any particular pipeline.

## The preflight

Three failures, each of which otherwise surfaces hours in, inside a fan-out task log:

1. **mhcmatch is not on the executor's `PATH`** — found at the first task rather than the twentieth.
2. **the installed release is not the one this overlay was written against.** The processes call CLI
   flags by name, so an older release exits 2 on an unknown flag deep in a task log. It has happened:
   PyPI served the previous release while the module had gained `--map-binder`, and both cassette
   tasks died with `unrecognized arguments`. Pin with `--mhcmatch_overlay_require_version`.
3. **the candidate table names no peptide or no allele column** — which, discovered later, reads as
   "this candidate named no allele we know", a real and *different* state.

In `off` mode it warns instead of failing, so an overlay that is merged but not yet enabled can
never break a run.

## Testing it

```bash
# topology, all five modes, no command run
nextflow run test/main.nf -stub-run -c overlay.config \
    --mhcmatch_overlay_mode cassette \
    --mhcmatch_cassette_unit_column context_peptide \
    --table ../../fixtures/S1.mhc1.candidates.tsv \
    --fasta ../../fixtures/S1.mhc1.windows.fasta
```

A stub runs no command, so it **cannot** catch an unrecognised CLI flag — exactly the class that has
bitten this integration before. Drop `-stub-run` (with mhcmatch installed) for the run that can.
