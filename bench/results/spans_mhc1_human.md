# Source-protein context discrimination, MHC-I (human + mouse)

For MHC-I **the peptide is the ligand** — there is no span to extend, so there is no span task. The
question instead is the antigen-processing one MHCflurry-2.0 poses (PMID 32711842): does an 8–11mer's
*source-protein context* look like a real ligand's?

```
python bench/bench_spans.py --spans spans.tsv \
    --proteome <p>/human.fasta.gz <p>/mouse.fasta.gz --cls mhc1 --seed 0
```

**Task.** AUROC separating 604,201 real eluted ligands **in their true source context** from
length-matched decoy windows drawn from the **same proteins**. Split by gene (475,317 train /
119,119 test; 24,294 / 6,098 genes); the model is fit on the train fold only.

## Result

| score | AUROC |
|---|---|
| full 12-position (flanks **+ the ligand's own termini**) | **0.814** |
| flank-only 6-position (upstream + downstream only) | 0.558 |
| shuffled-context control | 0.501 |

## Interpretation

**Read these two rows together — the gap between them is the whole point.**

The full 12-position score looks strong (0.814), but it is **not** a processing score. Six of its
twelve positions are the ligand's own first and last three residues, and for MHC-I those *are the
anchors* (P1–P3, PΩ). So most of that 0.814 is **binding**, re-measured. Quoting it as an
antigen-processing result would be a category error.

The honest processing signal is the **flank-only** score: **AUROC 0.558**. Real — the shuffled-context
control sits at 0.501, exactly chance — but weak. That is the expected answer, and it agrees with the
field: modern EL-trained binding predictors already absorb most of the processing signal, because
every eluted ligand's C-terminus *was* a proteasome product. What survives is a small residual, which
is exactly what MHCflurry-2.0's antigen-processing model extracts. `processing_score(flank_only=True)`
is that residual; use it to re-rank peptides the binding model already likes, never on its own.

## The class-I / class-II asymmetry falls out of the pooling check

Per-allele context PWMs, JSD against the pooled PWM:

| class | JSD (mean, top alleles by n) |
|---|---|
| MHC-II | **0.003 – 0.010** |
| MHC-I | 0.045 – 0.097 |

MHC-II context is essentially **allele-independent**, so pooling across alleles is justified — the
core is buried in the middle of the ligand and the termini are protease territory, not groove
territory. MHC-I context is **not** allele-independent, for exactly the reason above: the class-I
ligand's termini *are* its anchors, so a per-allele signal leaks into the ligand-internal positions.

This is a direct, quantitative confirmation of bind-first-trim-later, and it is why the class-I and
class-II entry points are separate functions (`processing_score` vs `presented_span`) rather than one
function with a `cls` flag.
