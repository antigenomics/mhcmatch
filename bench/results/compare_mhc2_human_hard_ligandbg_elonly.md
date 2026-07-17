# mhcmatch vs NetMHCIIpan-4.3i — human MHC-II, **eluted-ligand positives**

The same allele-specificity task as `compare_mhc2_human_hard_ligandbg.md`, with positives restricted
to (peptide, allele) pairs backed by at least one mass-spectrometry assay (`--el-only`, see
`bench/compare/provenance.py`). It answers *"can it find eluted ligands"* rather than *"can it
reproduce IEDB"*; the panel is EL-dominated but **not** EL-only, and the non-MS share is confounded
with allele.

**This is an evaluation stratum, not a training filter.** The model is fit on the whole corpus —
eluted ligands, binding assays, everything — exactly as it ships; `--el-only` changes only which
pairs are eligible to be *positives*. Training on the full corpus and tuning per task by parameter is
the house rule (`CLAUDE.md`), and binding-assay peptides are valid motif evidence: they do bind.

```
python bench/compare/run_compare.py --cls mhc2 --species human --benchmark holdout \
    --decoy-mode hard --background ligand --footprint adaptive --el-only
```

| stratum | metric | n alleles | mhcmatch | NetMHCIIpan-4.3i | Δ (mm−net) | p |
|---|---|---|---|---|---|---|
| medium | AUROC | 5 | 0.906 | **0.977** | -0.071 | 0.000 |
| medium | AUPRC | 5 | 0.639 | **0.730** | -0.091 | 0.000 |
| medium | PPV@P | 5 | 0.593 | **0.694** | -0.101 | 0.000 |
| frequent | AUROC | 20 | 0.901 | **0.952** | -0.050 | 0.000 |
| frequent | AUPRC | 20 | 0.558 | **0.682** | -0.124 | 0.000 |
| frequent | PPV@P | 20 | 0.521 | **0.666** | -0.145 | 0.000 |

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

| stratum | metric | Δ, all positives | Δ, EL positives |
|---|---|---|---|
| frequent | AUROC | -0.053 | **-0.050** |
| frequent | AUPRC | -0.124 | **-0.124** |
| medium | AUROC | -0.017 | **-0.071** |
| medium | AUPRC | -0.025 | **-0.091** |

Both tools score *higher* on eluted-ligand positives (mhcmatch frequent AUROC 0.892 → 0.901;
NetMHCIIpan 0.945 → 0.952) — an eluted ligand is a cleaner label than "bound in a competition assay"
— but the **gap barely moves** on frequent and *widens* on medium, where the surviving alleles are the
well-studied ones NetMHCIIpan knows best.

So the provenance stratum changes **what the number is about**, not who wins. The frequent-allele gap
survives it and remains the real one — that is where the work goes.

(A planning-stage estimate put the EL swing at +0.051 AUROC; that came from an ad-hoc harness with a
different split and decoy set. Measured properly here it is +0.009 for mhcmatch on frequent. The
smaller number is the right one.)

## Source-conditioning was tested and is not needed

The natural next move is an *adjusted general model* per provenance — one corpus, a `source`
parameter. The measured lever would be the core-offset prior, since EL boundaries are biological
(H/Hmax 0.720, peaked) while binding-assay boundaries are experimenter-chosen (0.990, flat as random
peptides), which suggests BA/in-silico queries should get a uniform prior. Held out, scoring each
query type with the corpus-learned prior vs a uniform one:

| query source | learned prior | uniform prior | Δ |
|---|---|---|---|
| EL | **0.922** | 0.912 | +0.010 |
| BA | **0.798** | 0.796 | +0.001 |

The learned prior helps eluted-ligand queries and is **harmless** on binding-assay ones, so a source
switch buys nothing (-0.001 at best). The general model already serves all three sources; the
existing tunables (`background`, `footprint`, `register`, `h`, `tau`) remain the per-task knobs. If
provenance ever enters the pmhc schema, re-test — but do not build the plumbing on spec.

## Why the join is on `(epitope, PMID)`

The pmhc schema carries no assay type, so provenance must come from the raw IEDB dump. Joining on
`(epitope, reference_id)` — present in both tables — avoids parsing the dump's single restriction
string (`HLA-DPA1*01:03/DPB1*04:01`) back into the alpha/beta pair `class2_key` consumes. It asks "in
this paper, was this peptide detected by mass spec?". A paper reporting MS for one allele and a
binding assay for another marks both EL; that is rare and errs toward keeping data. See
`bench/compare/provenance.py`.
