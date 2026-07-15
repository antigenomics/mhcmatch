# mhcmatch vs NetMHCpan-4.2b concordance — TESLA1 mhc1

Sample **TESLA1** (public), class mhc1. Common axis: presentation %rank (lower = stronger); Spearman ρ computed on a strength orientation so **ρ > 0 = agreement**. tier=full, background=proteome, footprint=adaptive, k-mer lengths=8,9,10,11, seed=0.

## Allele coverage

- scored by **both** (6): HLA-A02:01, HLA-A68:01, HLA-B15:07, HLA-B44:02, HLA-C03:03, HLA-C07:04
- mhcmatch only (0): —
- NetMHCpan-4.2b only (0): —
- neither (0): —

## View A — dense mhcmatch vs NetMHCpan-4.2b (all tiled k-mers × alleles)

- pooled Spearman ρ = **+0.726** over 41,754 (k-mer, allele) pairs
- strong-binder overlap (Jaccard): %rank≤0.5 = +0.247 (127 mm / 206 NetMHCpan-4.2b); %rank≤2.0 = +0.292 (528 mm / 883 NetMHCpan-4.2b)
- best-allele agreement (k-mers where either tool binds): **82%** (733/889)

### Per-allele Spearman ρ

| allele | ρ | n k-mers |
|---|---|---|
| HLA-A02:01 | +0.816 | 6,959 |
| HLA-B15:07 | +0.786 | 6,959 |
| HLA-B44:02 | +0.758 | 6,959 |
| HLA-C03:03 | +0.728 | 6,959 |
| HLA-A68:01 | +0.704 | 6,959 |
| HLA-C07:04 | +0.626 | 6,959 |

### Band confusion (rows mhcmatch, cols NetMHCpan-4.2b)

| mhcmatch＼NetMHCpan-4.2b | strong | weak | non-binder |
|---|---|---|---|
| strong | 66 | 33 | 28 |
| weak | 69 | 151 | 181 |
| non-binder | 71 | 493 | 40662 |

## View B — 3-way on the pipeline's own calls

Over 502 pipeline-called (epitope, best_allele) rows scored by all three:

| tool pair | Spearman ρ |
|---|---|
| mhcmatch~netmhc | +0.686 |
| mhcmatch~pipeline | +0.644 |
| netmhc~pipeline | +0.745 |

