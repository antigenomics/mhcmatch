# mhcmatch vs NetMHCIIpan-4.3i concordance — TESLA1 mhc2

Sample **TESLA1** (public), class mhc2. Common axis: presentation %rank (lower = stronger); Spearman ρ computed on a strength orientation so **ρ > 0 = agreement**. tier=full, background=proteome, footprint=adaptive, k-mer lengths=15, seed=0.

## Allele coverage

- scored by **both** (8): DRB1_1101, DRB1_1301, DRB3_0101, DRB3_0202, DRB4_0101, HLA-DPA10103-DPB10401, HLA-DQA10103-DQB10603, HLA-DQA10501-DQB10301
- mhcmatch only (0): —
- NetMHCIIpan-4.3i only (0): —
- neither (0): —

## View A — dense mhcmatch vs NetMHCIIpan-4.3i (all tiled k-mers × alleles)

- pooled Spearman ρ = **+0.559** over 16,072 (k-mer, allele) pairs
- strong-binder overlap (Jaccard): %rank≤0.5 = +0.170 (64 mm / 60 NetMHCIIpan-4.3i); %rank≤2.0 = +0.266 (224 mm / 247 NetMHCIIpan-4.3i)
- best-allele agreement (k-mers where either tool binds): **74%** (248/337)

### Per-allele Spearman ρ

| allele | ρ | n k-mers |
|---|---|---|
| DRB1_1101 | +0.717 | 2,009 |
| DRB1_1301 | +0.694 | 2,009 |
| HLA-DPA10103-DPB10401 | +0.681 | 2,009 |
| DRB3_0101 | +0.611 | 2,009 |
| DRB3_0202 | +0.590 | 2,009 |
| DRB4_0101 | +0.563 | 2,009 |
| HLA-DQA10501-DQB10301 | +0.341 | 2,009 |
| HLA-DQA10103-DQB10603 | +0.253 | 2,009 |

### Band confusion (rows mhcmatch, cols NetMHCIIpan-4.3i)

| mhcmatch＼NetMHCIIpan-4.3i | strong | weak | non-binder |
|---|---|---|---|
| strong | 18 | 16 | 30 |
| weak | 22 | 43 | 95 |
| non-binder | 20 | 128 | 15700 |

## View B — 3-way on the pipeline's own calls

Over 472 pipeline-called (epitope, best_allele) rows scored by all three:

| tool pair | Spearman ρ |
|---|---|
| mhcmatch~netmhc | +0.503 |
| mhcmatch~pipeline | -0.034 |
| netmhc~pipeline | -0.092 |

