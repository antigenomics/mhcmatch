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

---

# Which peptide should you actually use? (added after review)

The span model's point estimate is **not accurate enough to pick a peptide from**. Two independent
measurements say so, and they agree on the fallback: **use a fixed flank**.

## 1. Real crystals resolve a 13mer, not an 11mer

`bench/pdb_flanks.py` over the 93 pMHC-II structures of the Canonical2026 set (tcren-ms), assigning
the 9-mer core by groove burial:

| quantity | value |
|---|---|
| median **resolved** peptide length | **13** |
| median N-flank / C-flank beyond the core | **2 / 2** |
| structures resolving **≤11** residues | **13%** |
| structures resolving ≥2 flanking residues on both sides | 47% |

Length histogram: 9:1, 10:1, 11:10, 12:17, **13:25**, 14:11, 15:12, 16:11, 17:2, 19:2, 20:1.
(MHC-I control, same script: 180/260 are 9mers, 62 are 10mers — the method is sound.)

So **core ± 1 (11mer) is too short.** TCRmodel2 (PMID 37140040) and the fine-tuned AlphaFold of
Motmaen et al. (PMID 36802421) truncate their *input* to core ± 1, but that is a pipeline convention,
not a claim about what is ordered — 87% of real structures resolve more than 11 residues.
**Use core ± 2 (13mer)** for structure prediction: `ligand.STRUCTURE_FLANK`.

## 2. For an assay, coverage beats precision

Held-out (n=3000 cores), against the observed nested set. "Contains a full observed ligand" = the
emitted peptide brackets a real eluted span, i.e. an APC could produce the natural ligand from it.

| choice | both bounds ±1 | ±2 | ±3 | **contains a full observed ligand** |
|---|---|---|---|---|
| flank model (argmax span) | **31%** | 47% | 62% | 36% |
| fixed centred 13mer (core±2) | 26% | **51%** | 66% | 11% |
| fixed centred 15mer (core±3) | 28% | 50% | **79%** | 31% |
| fixed centred 17mer (core±4) | 19% | 39% | 62% | 52% |
| fixed centred 19mer (core±5) | 10% | 23% | 43% | 67% |
| fixed centred **21mer** (core±6) | 4% | 12% | 24% | **80%** |

**The model barely beats centring a 15mer on boundary error, and loses to it at ±3.** Its edge is the
exact-span hit rate (0.158 vs 0.069) — the question *"what was eluted?"*, not *"what should I make?"*.

For a CD4 assay the APC re-trims whatever you give it, so what matters is that the peptide **contains**
the natural ligand: a 21mer does so 80% of the time, the conventional 15mer only 31%. Longer also
tracks the MHC-II affinity optimum of ~18–20 aa (O'Brien 2008, PMID 19036163).
**Use core ± 6 (21mer)** for synthesis: `ligand.ASSAY_FLANK`.

## 3. Length calibration: the mean is right, the spread is too narrow

| | human | mouse |
|---|---|---|
| observed mean ligand length | 15.7 | 15.7 |
| predicted mean length | **15.5** | **15.6** |
| observed % at 15–16 | 33% | 38% |
| predicted % at 15–16 | **66%** | **66%** |
| observed % ≥18 | 22% | 19% |
| predicted % ≥18 | **4%** | **4%** |

The mean is well calibrated in both species, but the argmax **collapses onto the mode**: it emits a
15/16mer two-thirds of the time and almost never emits the long (≥18) or short (≤13) ligands that make
up ~39% of the real distribution. That is inherent to taking a point estimate from a broad
distribution — and it is the third reason not to treat the emitted span as *the* answer.

## Bottom line

| purpose | use | why |
|---|---|---|
| **structure** (TCRmodel2 / AF / Boltz) | **core ± 2 → 13mer** | median resolved crystal; 87% of structures resolve >11 residues |
| **synthesis / CD4 assay** | **core ± 6 → 21mer** | contains the true ligand 80% of the time; APC re-trims; affinity optimum 18–20 aa |
| **"what was actually eluted?"** | `presented_span(mode="modeled")` | 2.3× the exact-span hit rate of any fixed choice |
