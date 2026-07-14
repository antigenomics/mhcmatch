# Ligand-span recovery, MHC-II (human + mouse)

Given the binding core of a **held-out** eluted ligand, recover the ligand's observed span in its
source protein.

```
python bench/train_spans.py --iedb <iedb>/mhc_ligand_full.tsv.gz \
    --proteome <p>/human.fasta.gz <p>/mouse.fasta.gz --cls both \
    --out src/mhcmatch/data/ligand_context.tsv --spans-out spans.tsv
python bench/bench_spans.py --spans spans.tsv \
    --proteome <p>/human.fasta.gz <p>/mouse.fasta.gz --cls mhc2 --n-test 3000 --seed 0
```

**Data.** IEDB `mhc_ligand_full`, mass-spectrometry (eluted-ligand) assays only — binding-affinity
peptides have experimenter-chosen boundaries, which is exactly the label we are trying to predict.
5,571,576 assay rows → 1,918,444 class II → 1,323,569 EL → 373,904 spans after re-deriving every
coordinate (below). Split **by gene** (`GN=`), 295,869 train / 75,245 test, 14,295 / 3,582 genes.

**The model.** 12 terminus-relative context positions (3 upstream + the ligand's own first 3 + its
last 3 + 3 downstream — the NetMHCIIpan `-context` window) scored as log-odds against an order-1
Markov proteome null, plus a ligand-length prior. Allele-agnostic. **No free parameters**: the score
is `log P(L) + context log-odds`. (A tuned weight on the length prior looked better on the train
fold and did not transfer — 0.155 vs 0.158 held-out — so it was dropped, not shipped.)

## Result

| model | set-recall | IoU | median ΔN | median ΔC | mean ΔN | mean ΔC |
|---|---|---|---|---|---|---|
| modal length (15mer), centered on core | 0.069 | **0.792** | 2 | **1** | **1.81** | **1.93** |
| context only (no length prior) | 0.098 | 0.723 | 2 | 2 | 3.17 | 2.69 |
| **flank model (length + context)** | **0.158** | 0.765 | 2 | 2 | 2.25 | 2.31 |
| leak canary (proteins shuffled) | 0.048 | 0.730 | 2 | 2 | 2.58 | 2.62 |

*exact-span oracle ceiling **0.547*** — a core has **2.65 observed spans on average** (nested sets),
so for 45% of cores several spans are all legitimately correct and no model can pick "the" one.
Read 0.158 against 0.547, not against 1.0.

**set-recall** = top-1 span ∈ the observed nested set for that core. **ΔN/ΔC** = boundary error to
the nearest member of that set, **reported separately** — the whole hypothesis is N/C asymmetry, so
averaging the two sides would hide it.

## Interpretation

**The context signal is real, and both terms are needed.** The flank model more than doubles
exact-span recovery over centering a modal-length 15-mer (0.158 vs 0.069) — 29% of the achievable
ceiling. Length alone gives 0.069 and context alone 0.098; together 0.158, so they are complementary
rather than redundant. The **leak canary** (same model, protein sequences shuffled) collapses to
0.048, below the trivial baseline, confirming the gain comes from real source-protein context and
not from the length prior leaking through.

**But the model does not beat the trivial baseline on boundary error** (mean ΔN 2.25 vs 1.81, IoU
0.765 vs 0.792). It nails the exact observed span far more often, and is worse in the tail when it
misses. Centering is a strong prior — real ligands do sit roughly centred on their core — and a
model that is willing to move the boundaries pays for it on the occasions it moves them wrongly.
Reported rather than buried: if you need a robust default and not the single most likely ligand,
centring a 15-mer is still a defensible choice.

This modesty is expected. NetMHCIIpan reports that its `-context` encoding adds only ~1–2 AUC0.1
points on eluted-ligand data (PMID 30446001), and that it **degrades** performance on CD4 T-cell
epitope benchmarks (PMID 32406916). Flank models predict *ligands*, not *immunogenicity*.

## Two corrections this benchmark forced

1. **IEDB's annotated coordinates cannot be used as labels.** They are wrong for ~8.8% of rows and
   *silently* wrong for ~3.8% — the peptide substring-matches its protein while the annotated start
   points somewhere else entirely (signal-peptide / isoform numbering). Every coordinate here is
   **re-derived** by exact substring match against the same FASTA used at inference, keeping only
   unique occurrences (2,925 zero-occurrence and 700 multi-occurrence spans dropped).
2. **A gene split is not sufficient on its own.** The same peptide occurs in several genes
   (paralogs, shared domains), so 2,790 train spans shared a peptide with a test gene. The
   split-integrity assertion caught it; they are purged from train.

## Known-biology control

Proline is the aminopeptidase stop signal — **enriched inside the ligand, depleted in the flank**:

| position | Pro vs proteome |
|---|---|
| `ligN+2` (inside the ligand) | **2.00×** |
| `flankN-1` (in the flank) | **0.25×** |

Note the sign: a model that put a proline *preference* in the flank would have it backwards. Half
the signal sits inside the peptide — which is precisely why the context window is 3+3+3+3 and not
flanks alone. Asserted in `tests/test_mhcmatch.py::test_vendored_span_table_recovers_known_biology`.

**Allele pooling is justified, not assumed.** Per-allele context PWMs are near-identical to the
pooled one (JSD 0.003–0.010 across the 7 class-II alleles with ≥5k ligands) — trimming is protease
biology, not groove biology. Pooling also unlocks the ~70% of class-II EL records whose restriction
is only a placeholder (`HLA class II`).

**Cysteine is clamped.** C is depleted 8–11× at the ligand termini but *not* in the flanks — the
flanks are not in the detected peptide, so this is mass-spectrometry chemistry, not processing. Left
in, the model would refuse every Cys-containing ligand.

## Positive control (public)

**MBP85-99 / DRB1\*15:01.** From core `VHFFKNIVT` in human MBP, the model returns
`NPVVHFFKNIVTPR` — the canonical DR2 ligand `ENPVVHFFKNIVTPR` minus one N-terminal residue
(ΔN = 1, ΔC = 0), with the true span among the scored alternatives.
