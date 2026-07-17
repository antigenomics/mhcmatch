# mhcmatch vs NetMHCIIpan-4.3i — human MHC-II, **eluted ligands only**

The same allele-specificity task as `compare_mhc2_human_hard_ligandbg.md`, restricted to
(peptide, allele) pairs with at least one mass-spectrometry assay (`--el-only`, see
`bench/compare/provenance.py`). This is the number a *presentation* predictor should be judged on:
the panel is EL-dominated but **not** EL-only, and the non-MS share is confounded with allele.

```
python bench/compare/run_compare.py --cls mhc2 --species human --benchmark holdout \
    --decoy-mode hard --background ligand --footprint adaptive --el-only
```

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | 95% CI | p |
|---|---|---|---|---|---|---|---|
| medium | AUROC | 5 | 0.930 | **0.977** | -0.047 | [-0.074, -0.020] | 0.001 |
| medium | AUPRC | 5 | 0.638 | **0.751** | -0.113 | [-0.203, -0.024] | 0.000 |
| medium | PPV@P | 5 | 0.564 | **0.666** | -0.102 | [-0.221, -0.020] | 0.001 |
| frequent | AUROC | 20 | 0.904 | **0.951** | -0.046 | [-0.083, -0.014] | 0.000 |
| frequent | AUPRC | 20 | 0.564 | **0.676** | -0.112 | [-0.183, -0.049] | 0.000 |
| frequent | PPV@P | 20 | 0.534 | **0.650** | -0.116 | [-0.186, -0.059] | 0.000 |

## The headline: there is no rare stratum

**23 of 52 human MHC-II panel alleles cannot support this benchmark at all.** Fifteen have **zero**
eluted ligands — DRB1\*01:02, DRB1\*01:03, DRB1\*03:02, DRB1\*03:03, DRB1\*04:03, DRB1\*08:02,
DRB1\*09:01, DRB1\*11:02, DRB1\*11:03, DRB1\*11:04, DRB1\*13:03, DRB1\*16:02, DRB4\*01:01,
HLA-DQA1\*03:01-DQB1\*03:01, HLA-DQA1\*05:01-DQB1\*03:02 — and eight more fall under a 20-ligand
floor (DRB1\*04:02 has 3, DRB1\*04:04 has 2, DRB1\*13:01 has 1).

The `rare` stratum disappears entirely. **So the "mhcmatch wins rare AUROC" result in
`compare_mhc2_human_hard_ligandbg.md` (and in `ROADMAP.md` §6) is measured on alleles that have no
eluted ligands** — it is a binding-assay benchmark reported as a presentation one. That needs saying
in the paper, and it is the single most load-bearing thing in this file.

## Effect on the gap: essentially none

| stratum | metric | all-provenance Δ | EL-only Δ |
|---|---|---|---|
| frequent | AUROC | -0.052 | **-0.046** |
| frequent | AUPRC | -0.125 | **-0.112** |
| medium | AUROC | -0.016 | **-0.047** |
| medium | AUPRC | -0.025 | **-0.113** |

Both tools score *higher* on EL-only positives (mhcmatch frequent AUROC 0.893 → 0.904; NetMHCIIpan
0.945 → 0.951) — eluted ligands are simply a cleaner label than "bound in a competition assay" — but
the **gap barely moves** on frequent (-0.052 → -0.046) and *widens* on medium, where the surviving
alleles are the well-studied ones NetMHCIIpan knows best.

So provenance filtering is a **correctness fix, not a gap-closer**, exactly as measured during
planning: it changes which alleles are honestly evaluable and what the number means, not who wins.
The frequent-allele gap survives it and remains the real one.

(A planning-stage estimate put the EL-only swing at +0.051 AUROC; that came from an ad-hoc harness
with a different split and decoy set. Measured properly here it is +0.011 for mhcmatch on frequent.
The smaller number is the right one.)

## Why the join is on `(epitope, PMID)`

The pmhc schema carries no assay type, so provenance must come from the raw IEDB dump. Joining on
`(epitope, reference_id)` — present in both tables — avoids parsing the dump's single restriction
string (`HLA-DPA1*01:03/DPB1*04:01`) back into the alpha/beta pair `class2_key` consumes. It asks "in
this paper, was this peptide detected by mass spec?". A paper reporting MS for one allele and a
binding assay for another marks both EL; that is rare and errs toward keeping data. See
`bench/compare/provenance.py`.
