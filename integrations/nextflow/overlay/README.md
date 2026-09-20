# mhcmatch as an overlay on a neoantigen pipeline you already run

For a host pipeline that already does variant calling, HLA typing and expression quantification, and
reaches a candidate table with those in hand — the **nf-core/sarek → VEP → pVACtools** shape, with
**OptiType** and **arcasHLA / HLA-LA** for the typing. It attaches at named seams and **contributes
zero new scoring processes**: everything is an alias of `../mhcmatch/`.

## Why an overlay rather than copying the processes in

Because a copy is a fork, and a fork drifts. A hand-wired integration we reviewed had drifted five
ways within a few releases:

| drift | consequence |
|---|---|
| its `-stub-run` hard-typed an 18-column header | stubs no longer matched the real 57-column table — the exact drift our own module repaired in 2026-09-20 and now prevents by asking the library for its schema |
| it emitted `native_table` where ours emits `native_tsv` | a channel rename that only shows up when a downstream `include` is added |
| its image was pinned to mhcmatch **1.0.1** | sixteen releases behind, and nothing in a run says so |
| its `rank_threshold` defaulted to **2.0** | the weak class-I cut, applied at the predictor — candidates its own later aggregation could still have ranked are gone before they reach it. The same number is the ***strong*** cut for class II: measured on one window pair against `DRB1*15:01`, **0 of 56** scored pairs survive it. That fork passed 100 on its class-II path, so the class-II hazard was latent rather than active — one `include` away |
| its class-II bridge called `mhcmatch predict` **once per allele**, in a Python loop | `predict` takes the whole allotype list and batches internally. N sequential calls also reload the model and rebuild the calibration N times |

An alias cannot drift. A sixth divergence will.

## Wiring it in

**Level 0 needs no edit to your repository at all.** One `-c`, plus your own predictor switch if you
have one — a config file cannot reach inside your `if`:

```bash
nextflow run . -c /path/to/mhcmatch/integrations/nextflow/overlay/overlay.config \
    <your predictor switch> mhcmatch
```

It pins the image every mhcmatch-invoking process runs in and sets `--rank-threshold none`, so
filtering happens where you can see what was filtered. It changes nothing else.

**A `withName` selector matches the alias.** If you `include { MHCMATCH_PREDICT as MHCMATCH_MHCI }`
— the ordinary way to serve two classes from one process — then a selector naming `MHCMATCH_PREDICT`
binds to nothing and says nothing. So the image pin uses a deliberately broad literal,
`.*MHCMATCH_[A-Z0-9_]*$`: it sets `container` and `conda` and nothing else, and pinning the image on
a process that never calls mhcmatch costs nothing. To bind other names, add your own `withName`
block — process directives compose. It cannot be a parameter: Nextflow rejects an interpolated
`withName:` and the whole config then fails to parse.

**Levels 1–3** need two lines. Nextflow has no way to interpose a process on a channel from config,
so this part is a script edit — but it is +1 include and a 7-line wrap.

In your config — **one line**, because `overlay.config` pulls in the module's own:

```groovy
includeConfig '/path/to/mhcmatch/integrations/nextflow/overlay/overlay.config'
```

In your neoantigen workflow, beside your other includes:

```groovy
include { MHCMATCH_OVERLAY } from '/path/to/mhcmatch/integrations/nextflow/overlay/overlay.nf'
```

then **one call, three arguments**, placed ahead of your own selector:

```groovy
MHCMATCH_OVERLAY(
    SCORED_CANDIDATES.out.csv,          // [ meta, table, cls ] -- your ranked candidate table
    MERGED_WINDOWS.out.sequences,       // [ meta, fasta, cls ] -- the mutation windows
    ch_alleles                          // [ meta, 'HLA-A*02:01,...' ] or [mhc1:'…', mhc2:'…']
)
versions = versions.mix( MHCMATCH_OVERLAY.out.versions )
```

and at the assignment of the FASTA your construct optimiser consumes — a **ternary, not a cut**:

```groovy
vaccine_aminoacid = params.mhcmatch_overlay_mode == 'cassette'
                        ? MHCMATCH_OVERLAY.out.vaccine
                        : YOUR_SELECTOR.out.concatemer
```

**There is deliberately no fourth argument.** An earlier signature took the FASTA your construct
step emits, purely to hand it back in the four modes that do not replace it — which is a channel
cycle the moment you use both seams: the overlay would need your selector's *output* while producing
its *input*, and Nextflow cannot schedule that. `out.vaccine` is empty except in `cassette` mode, so
the ternary above keeps your own construct, and it is the same shape as any construct-variant
switch you already have.

**Your own construct step still runs**, in every mode including `cassette`. Cutting it would take
its other outputs down with it — an epitope-marking step, a neoantigen table, a dashboard — which is
a silent degradation of a deliverable and precisely what this must not do.

### One thing to set before `cassette` mode

**The cassette unit must be the long (~27 aa) window, and on this path there is no window FASTA to
rebuild it from**, so name the column of your candidate table that carries it:

```bash
--mhcmatch_vector_unit_column context_peptide     # whatever your table calls it
```

`MHCMATCH_CASSETTE` **refuses to start** without either that or a `--context` FASTA, rather than
falling back to `peptide` — which on a scored table is the *minimal* epitope, and a 9-mer loads onto
any cell without costimulation: the tolerising configuration (PMID 17911588). It used to default to
one upstream's spelling, which answered for that upstream and silently built the tolerising cassette
for every other.

### The second seam, and it is the one `rerank` needs

**Your selector sorts on a column it names, so publishing ours beside it changes nothing.** A
selector picks the peptides that enter a construct by sorting on a column whose name it fixes —
typically a literal `score` — and `rank pairs --passthrough --prefix mm_` adds `mm_score` and leaves
that column exactly as it found it. That is the right default and it is also why `annotate` is a
genuine no-op on every deliverable. Set `params.mhcmatch_overlay_score_column` if yours is named
something else.

`rerank` mode writes `score`, and it goes back to **your** selector — your linkers, your length
budget, your allele groups, our ordering. Feed the second emitted channel where you fed your own
scored table:

```groovy
ch_scored = MHCMATCH_OVERLAY.out.scored          // [ meta, table, cls ]
YOUR_SELECTOR(
    ch_scored.filter  { m, t, cls -> cls == 'mhc1' }.map { m, t, cls -> [ m, t ] }
    .join( ch_scored.filter { m, t, cls -> cls == 'mhc2' }.map { m, t, cls -> [ m, t ] } )
    …                                            // the rest of your own join, unchanged
)
```

In every mode but `rerank` that channel is the object you passed in, verbatim — so the wrap is safe
to merge before anyone has decided anything, exactly like the ternary.

**Every multiplicative correction you applied is preserved, because they are clinical and this model
does not contain them.** A scorer of this shape builds `score = base × driver × recurrence ×
hla_freq × hla_loss`, and only `base` is a question about peptides. A driver boost of **100.0** and
an HLA-loss weight of **0.0** for an allele the tumour has lost are both ordinary settings — the
first deliberately two orders of magnitude, the second a veto. So `rerank` replaces the base and
carries the factors:

```
score := W × minmax(mm_score) × driver(driver_class) × hla_weight × recurrence_weight
         × hla_frequency_weight
```

`W` is the sum of your own four base weights, read from your params when they are loaded in the same
session, so the numbers stay in the range your operators read. A factor whose column is absent is
1.0 — a run with no recurrence library has no `recurrence_weight` and is scored as your own scorer
would have scored it. Set `params.mhcmatch_overlay_keep_weights` to your own factor columns. Your
value is kept as `<column>_host`, `score_host` by default; nothing is destroyed.

Measured on the 502 class-I candidates of the public TESLA1 fixture, the two orderings agree at
Spearman **ρ = 0.666**, sharing **4 of the top 20** and **15 of the top 50** — so this changes the
construct, which is the point of the mode.

## The five modes

| `--mhcmatch_overlay_mode` | what runs | what it touches |
|---|---|---|
| **`off`** *(default)* | preflight only | nothing. `scored` is your object verbatim and `vaccine` is empty, so both wrappings are inert. Safe to merge before anyone has decided anything |
| `annotate` | the whole rerank arm — reranked table, neoantigen table, units, cassette and cohort score — published **beside** your table under `${outdir}/mhcmatch/` | **nothing of yours.** Every deliverable you own is byte-identical; you get our columns next to yours for comparison. The mode to start in. It is additive, not cheap: the arm runs the cassette build and its safety screen |
| `rerank` | as `annotate`, then `score` is rewritten from this model and handed back | **your selector, our ordering.** Your linkers, length budget and allele groups all still apply |
| `denovo` | our de novo arm on your merged window FASTAs | a third table; orthogonal to your own predictor choice |
| `cassette` | rerank → select → assemble → the cohort score | our selector replaces yours; the mode that changes what your optimiser receives |

`annotate` and `rerank` differ by exactly one process, and that process is the point of the mode.

Deliberately **not** a new value of your own construct-schema parameter: that would be a second diff
in a second place (your validator), and it conflates *what the construct contains* with *who designed
it*. Those are orthogonal and both are worth varying.

## Failsafe

Four mechanisms, in the order they fire.

**1. `MHCMATCH_PREFLIGHT` runs before anything, in every mode.** It computes nothing. It checks that
`mhcmatch` is on PATH, that the installed release matches `--mhcmatch_overlay_require_version`
(**1.20.0**), and that the candidate table names a peptide column and an allele column — reading the
header with the delimiter it actually has, comma or tab, because a candidate table is either.
Each failure names its cause.

Without it, each of those surfaces hours later inside a fan-out task log: a missing install as
`command not found` at the twentieth task, a version mismatch as `unrecognized arguments` (which has
happened — PyPI served the previous release while the module had gained `--map-binder`, and both
cassette tasks died), and a missing column as an empty `allele` field, which reads as *this candidate
named no allele we know* — a real and different state.

In `off` mode preflight **warns rather than fails**, so a merged-but-unused overlay can never break a
run.

**2. Additive output.** The overlay writes only under `${params.outdir}/mhcmatch/` and rebinds no
channel you read. Its `publishDir` selector is **enumerated rather than `.*MHCMATCH_.*`** — a
catch-all also matches *your* mhcmatch-named processes and relocates their published files, which is
not additive whatever the directory is called.

**2a. `cassette` mode refuses to start when a quota is set.** With `--mhcmatch_vector_quota` the
assembly step emits two records, `cassette_composed` and `cassette_topk`, so the comparison is on
your own candidates — and a codon optimiser handed two records brackets each one and writes two
constructs under one sample name. Run the comparison under `rerank` and read the FASTA yourself.

**3. `cassette` mode wraps rather than cuts**, as above.

**4. The error strategies are asymmetric on purpose.** Annotation-only processes are `ignore` — a
failed mimicry step should cost one missing extra table, not your cohort. The cassette path is
`terminate`, and the reason is not symmetric: under `ignore` a failed upstream process empties its
channel, and **every process downstream of an empty channel is skipped without error**. `ignore` on
the cassette path would produce a run that reports success and builds nothing.

Forgetting the module's params used to be its own trap, and a bad one: every `params.mhcmatch_*`
would be undefined, undefined is falsy, and `mhcmatch_vector_screen` therefore read as **off** — a
cassette built with no safety screen, reported as one WARN among a dozen. `overlay.config` now
includes the module's config itself, so the ordinary path cannot get it wrong; `overlay.nf` still
refuses to run when the params are missing, because a caller can always include this file the wrong
way.

## Verifying it, without a reference bundle and without patient data

**V1 — the stub, all five modes.** Seconds, no reference genome, no install needed:

```bash
for m in off annotate rerank denovo cassette; do
    nextflow run test/main.nf -stub-run -c overlay.config \
        --mhcmatch_overlay_mode $m \
        --table ../../fixtures/S1.mhc1.candidates.tsv \
        --fasta ../../fixtures/S1.mhc1.windows.fasta
done
```

It proves the channel topology resolves in every mode, that `scored` carries **the object that was
passed in** in four of five, and that `vaccine` carries something in exactly one.

**A stub runs no command, so V1 structurally cannot catch an unrecognised CLI flag** — the class that
has actually bitten this integration. It also cannot reach the guards that live in a `script:`
block, which is why `--mhcmatch_vector_n0` and `--mhcmatch_vector_unit_column` are not needed for
V1 and are needed for V2. That is what V2 is for.

**V2 — a real run on the fixture.** `../../fixtures/` holds two synthetic donors — a samplesheet,
tab-separated candidate tables carrying pVACseq's own spellings (`MT Epitope Seq`, `HLA Allele`,
`Gene Name`, `WT Epitope Seq`) plus the 27-mer window as `context_peptide`, window FASTAs with the
colon-delimited header, and an OptiType-shaped typing file — with invented sequences and **no real
variant or identifier of any kind**. That is deliberate: the fixture has to exercise the real shapes
and it must not be a patient's.

What a real run proves that the stub cannot:

1. every CLI flag the overlay passes is accepted by the installed release;
2. the table is read with the delimiter it actually has — a tab-only header split turns a wide CSV
   schema into one column, and six readers in this library had that bug;
3. all your columns survive `--passthrough` under their own names, in your order, ahead of ours.
   Measured on one 57-column comma-separated candidate table: **91 columns out, the caller's 57
   preserved in the caller's order**;
4. a column of yours that collides with one we emit is preserved as `<name>_in` rather than
   overwritten — measured on a table carrying both `score` and `group`, which both collide. The
   synthetic fixtures carry neither, so this one needs a table of your own;
5. `--mhcmatch_vector_unit_column context_peptide` resolves, so the cassette is built from the
   **27-mer window** and not from the minimal epitope. A 9-mer loads onto any cell without
   costimulation and is the tolerising configuration; neither arm is allowed to inject one;
6. under `rerank`, the rewritten table still names `score` and is readable by a host selector —
   57 of its own columns plus ours, not one column and a `KeyError`.

**V3 — the rescore transform, on its own.** It is the one step whose output another program parses
by column name, so it is checked against a real reranked table rather than only through the DAG:

```bash
mhcmatch rank pairs <your>.candidates.tsv --cls mhc1 --passthrough --prefix mm_ \
    --out /tmp/reranked.tsv
sed -n "/^import csv, sys$/,/^PY$/p" rescore.nf | sed '$d' > /tmp/rescore.py
python3 /tmp/rescore.py /tmp/reranked.tsv /tmp/out.csv mm_score 7.0 100.0 100.0 \
    hla_weight,recurrence_weight,hla_frequency_weight
```

On the TESLA1 fixture that is 502 rows in and 502 out, 92 columns, 0 unscored, and the ρ = 0.666
above.

**What this cannot verify here, stated rather than papered over:** a host pipeline's own upstream
legs. Those need a reference bundle measured in hundreds of gigabytes and real sequencing input, so
the overlay is verified at the module boundary and integration against a live host is the operator's
step.
