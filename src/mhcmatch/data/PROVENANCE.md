# Vendored data provenance

## `mhci_pseudo.fa` / `mhcii_pseudo.fa`

NetMHCpan-style **34-residue MHC pseudosequences** — the polymorphic groove positions that contact
the peptide (class I: α1/α2 of the MHC heavy chain; class II: α1 of the α-chain + β1 of the
β-chain). `X` marks an ambiguous/unresolved position.

Alleles sharing a 34-mer are collapsed to one record. The header lists **every** allele of the
group, space-separated: `>ALLELE [ALLELE ...]|n=<count>`. All of them are keys in
`pseudoseq.load_pseudo`, so a query for any allele in the group returns that group's sequence —
which *is* that allele's own sequence, since the group is defined by exact 34-mer identity.

- **MHC-I:** 5407 unique sequences over **20082 alleles** (human HLA-A/B/C/E/F/G, mouse H-2, others).
- **MHC-II:** 2209 unique sequences over **11048 alleles** (human HLA-DR/DQ/DP, mouse H-2 I-A/I-E, others).

Two sources, built by the sibling `antigenomics/tcren` repo's `scripts/build_pseudo_fasta.py`:

1. **NetMHCpan's tables** (`MHC_pseudo.dat`, `pseudosequence.2023.all.X.dat`) — 12997 MHC-I alleles.
   Authoritative wherever present.
2. **IPD-IMGT/HLA 3.65.0** (`ANHIG/IMGTHLA`, `alignments/{A,B,C,E,F,G}_prot.txt`) — **+7085**
   class-I alleles the NetMHCpan table has never covered. It lags IMGT and omits **HLA-F entirely**.
   The 34 groove positions are *not* hardcoded: they are recovered by consensus from the alleles the
   table already knows, cross-checked between genes (HLA-B and HLA-C solve independently and agree),
   then applied to genes with too few knowns by aligning reference sequences — with **HLA-E and
   HLA-G as positive controls** (both round-trip 100%, which is what licenses HLA-F, whose 0 known
   alleles leave nothing to check directly). Verified by re-deriving every known allele:
   **21935 exact, 4 mismatch (0.018%)**. The 4 are indel-bearing alleles (A\*24:164, A\*24:399,
   A\*32:80, B\*51:50) where NetMHCpan-4.2 places the gap one slot from IMGT 3.65.0; NetMHCpan wins
   every conflict, so no already-covered allele can change. 81% of the added alleles simply join an
   existing 34-mer group — new HLA alleles usually differ outside the groove.

Used by `mhcmatch.pseudoseq` as the allele-similarity alphabet for the cross-allele diffusion model
(see `../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` §4). Regenerate with:

    for g in A B C E F G; do
      curl -sSo ~/vcs/tmp/imgt/${g}_prot.txt \
        https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/${g}_prot.txt
    done
    python ../tcren-ms/scripts/build_pseudo_fasta.py \
        --mhci  ~/work/academy/software/netMHCpan-4.2/data/MHC_pseudo.dat \
        --mhcii ~/work/academy/software/netMHCIIpan-4.3/data/pseudosequence.2023.all.X.dat \
        --imgt-alignments ~/vcs/tmp/imgt \
        --out src/mhcmatch/data

These files are static reference data and small (~800 KB total), so they are vendored rather than
fetched. Re-sync from `tcren` if the pseudosequence definition is updated upstream.

**History (2026-07-16).** Until this date the header carried only the group's *first* allele, so the
other 8854 of MHC-I's 12997 alleles (68%) — and 8839 of MHC-II's 11048 (80%) — were **silently
unresolvable**, among them common specificities like HLA-B\*14:02, B\*18:05 and C\*03:04. The
collapse was always correct; only the name index was lost. Restoring it left every 34-mer
byte-identical (asserted at regeneration) and lifted the MixMHCpred3 benchmark from maxF1 0.8807 to
0.8908. Both this file and `tcren`'s builder were fixed; a re-sync from an unfixed `tcren` would
silently reintroduce the bug. The IMGT source was added at the same time, taking the human MHC-I
reference panel from 166/203 scorable alleles to **203/203**.

## `mhc2_alpha_prior.tsv`

`DP/DQ beta chain -> most likely alpha chain`, used by `pseudoseq.class2_key(..., impute_alpha=True)`
to key class-II records that type only the β chain. **4824 human MHC-II records (1.5% of the panel;
2516 of them HLA-DPB1\*11:01 alone) carry an empty `mhc_a`** — they arrive as the groove-less key
`-DPB11101` and were previously dropped by `Store.from_records` outright. DRA is monomorphic so DR
needs no table (it is hardcoded in `class2_key`); DPA1/DQA1 are polymorphic and must be learned.

**Derived** (not experimental) from the IEDB-derived pmhc panel — see `bench/compare/SOURCES.md` for
the panel itself. A β is listed only when its **34-mer groove** is ≥95% determined over ≥50
fully-typed ligands. 20 β chains qualify; 9 are refused.

The criterion is the *groove*, not the allele name nor its 2-digit group — measured, not assumed.
`DQA1*01:02` and `DQA1*01:05` share the group `DQA1*01` but have **different 34-mers**, so a
group-level rule reads 100% certain on `DQB1*05:02` while the sequence is still a 58/42 coin flip
(DQA1's polymorphism sits in the α1 domain the pseudosequence samples). Refused β chains keep their
α-less key and stay unscorable on purpose: a wrong groove scores silently, which is worse than not
scoring. The table rediscovers DQ2.5 (`DQB1*02:01`–`DQA1*05:01`) and DQ8 (`DQB1*03:02`–`DQA1*03:01`)
blind, which is the expected linkage disequilibrium and a check that the method is sane.

**Where it is on, and where it is not — both measured.** The *lookup* path
(`class2_from_name`, `class2_key`) imputes by **default**: querying `HLA-DPB1*11:01` returned `nan`
before and now resolves to a real groove, which is a strict win. The *panel* path
(`Store.from_records`/`from_pmhc`) defaults to **off**, i.e. beta-only records stay dropped as they
always were. Admitting them to the reference panel was tested over the 13 alleles whose reference set
grows and it does **not** help — held-out AUROC **−0.0019**, AUPRC **−0.0012**, and the damage scales
with the merge: `HLA-DPA10201-DPB11101` gains 2339 ligands (+89%) and loses **0.0155 AUROC**. A study
that skipped α-typing produced noisier ligand calls too, so a missing α marks data quality, not just
absent metadata. (Caveat: the held-out positives are fully-typed ligands, so this measures whether
orphan ligands predict *typed* ones. Whether they predict orphan-like ligands is untested.)

If turned on for the panel (`tier=full`): stranded ligands 4782 → 635, `HLA-DPA10201-DPB11101`
2618 → 4957. Regenerate the table:

    # see the snippet in this file's git history (commit that added mhc2_alpha_prior.tsv);
    # it groups panel ligands by beta and keeps the modal 34-mer where P>=0.95 and n>=50.

**Known gap.** 9 DQ β chains fail the bar and stay unresolvable (635 ligands), e.g. `-DQB10503`
(P(groove)=54%) and `-DQB10502` (n=12). Rare DQ haplotypes are not in tight enough LD to impute.

**History (2026-07-16).** Until this date the header carried only the group's *first* allele, so the
other 8854 of MHC-I's 12997 alleles (68%) — and 8839 of MHC-II's 11048 (80%) — were **silently
unresolvable**, among them common specificities like HLA-B\*14:02, B\*18:05 and C\*03:04. The
collapse was always correct; only the name index was lost. Restoring it left every 34-mer
byte-identical (asserted at regeneration) and lifted the MixMHCpred3 benchmark from maxF1 0.8807 to
0.8908. Both this file and `tcren`'s builder were fixed; a re-sync from an unfixed `tcren` would
silently reintroduce the bug.

## `ligand_context.tsv`

The ligand-span (flank/context) model consumed by `mhcmatch.ligand.load_span_model()`. Per class
(`mhc1`, `mhc2`): residue frequencies at the 12 terminus-relative context positions — 3 upstream in
the source protein, the ligand's own first 3 and last 3, and 3 downstream (the NetMHCIIpan
`-context` window, PMID 30446001) — plus a ligand-length prior. Laplace-smoothed at fit time, so the
runtime carries no smoothing parameter. Allele-agnostic (justified, not assumed: per-allele context
PWMs are within JSD 0.003–0.010 of the pooled one for MHC-II).

Fit from **IEDB** `mhc_ligand_full` (mass-spectrometry / eluted-ligand assays only — binding-affinity
peptides have experimenter-chosen boundaries, which is the very label being modelled) against the
UniProt reference proteomes **UP000005640** (human) and **UP000000589** (mouse). 604,201 MHC-I and
373,904 MHC-II spans survive; every coordinate is **re-derived** by unique exact substring match
rather than trusting IEDB's annotated `Starting Position`, which is wrong for ~8.8% of rows and
*silently* wrong for ~3.8%. Cysteine's log-odds at the ligand-terminal positions is clamped to the
proteome background: C is depleted 8–11× there but not in the flanks, i.e. it is mass-spectrometry
chemistry (alkylation / missed ID), not processing biology.

Inputs are distributed via the public HF dataset
<https://huggingface.co/datasets/isalgo/pmhc_data> (`dump/mhc_ligand_full.tsv.gz`,
`proteome/human.fasta.gz`, `proteome/mouse.fasta.gz`). Regenerate with:

    python bench/train_spans.py \
        --iedb <pmhc_data>/dump/mhc_ligand_full.tsv.gz \
        --proteome <pmhc_data>/proteome/human.fasta.gz <pmhc_data>/proteome/mouse.fasta.gz \
        --cls both --out src/mhcmatch/data/ligand_context.tsv

Held-out validation: `bench/results/spans_mhc2_human.md`, `bench/results/spans_mhc1_human.md`.

## `aa_tables.py`

Amino-acid property scales for `mhcmatch.immuno`: `DESCRIPTORS` (17 families, 102 components),
`HYDROPHOBICITY` (45 scales), `MJ_PARTITION`, and the two-component property basis
`PROPERTY_PC1` / `PROPERTY_PC2`. **Generated, never transcribed** — the vendoring script reads the
upstream tables and emits a literal Python module, so a re-run reproduces the file byte-for-byte and
a typo is impossible.

Two upstreams, both GPL-3.0-or-later like mhcmatch, so copying is licence-clean:

- **`peptides` 0.5.0** (<https://github.com/althonos/peptides.py>, PyPI), file
  `peptides/tables/__init__.py` → `DESCRIPTORS` + `HYDROPHOBICITY`. **Derived/computed**: each
  family is a published descriptor scale fitted by its own authors, redistributed by the package.
  It is *not* a dependency — nothing at runtime imports it.
- **`tcren` 2.8.0** (`antigenomics/tcren`), file `src/tcren/data/MJ1985_partition_energies.csv`
  = AAindex `MIYS850101` → `MJ_PARTITION`. **Derived/computed**: a statistical contact potential
  Miyazawa & Jernigan fitted to a PDB structure set. Larger = more hydrophobic.
  Re-fetch: <https://www.genome.jp/dbget-bin/www_bget?aaindex:MIYS850101>.

`PROPERTY_PC1` / `PROPERTY_PC2` have no upstream — they are **derived/computed** from the two tables
above: the first two principal components of the 20 × 142 property matrix (residues × scales),
column-standardized over the 20 residues, by SVD. **Label-free**, so identical under every
leave-one-dataset-out refit; sign convention is that each component's largest-magnitude residue score
is positive. PC1 carries **32.79 %** of total variance and is a hydropathy axis with residue order
`I F L W V M C Y A P G T H S Q N E K D R`; PC1 + PC2 carry **51.2 %**, and 10 components carry
**91.3 %** (`bench/results/ipred_pca.md`). Regenerate with `python bench/ipred/pca.py`. They were
first vendored inside a retired artifact as `residue_scores`; it left the wheel and these vectors,
which never depended on its fit, stayed.

Regenerate (from the `2026-mhcmatch-benchmark` repo; `peptides` must be downloaded and unpacked
first — it is not installed in any venv here):

    pip download peptides==0.5.0 --no-deps -d /tmp/pep && unzip /tmp/pep/peptides-0.5.0*.whl -d /tmp/pep/pkg
    python bench/immuno/vendor_aa_tables.py \
        --peptides /tmp/pep/pkg --tcren ~/vcs/code/tcren \
        --out src/mhcmatch/data/aa_tables.py

Verified at generation: 17 families / 102 components / 45 hydrophobicity scales / 20 MJ residues;
every table covers exactly AA20; Kidera KF1(A) = −1.56, VHSE1(A) = 0.15, Kyte-Doolittle I/R =
4.5/−4.5, MJ A/F = 2.36/4.37 all match the published values. Non-standard residues (B, J, O, U, X, Z)
are absent by construction so a caller must decide how to handle them rather than silently scoring 0.

## `contact_profile.py`

Per-position **TCR↔peptide** contact frequency by (MHC class, peptide length), backing
`mhcmatch.immuno`'s `"contact"` anchor scheme — a continuous positional weighting that needs no
anchor call at all. Generated, not hand-written.

Built from three `tcren` 2.8.0 files (`antigenomics/tcren`):

| input | path in `tcren` | what is taken | provenance |
|---|---|---|---|
| `contacts_2026.csv` | `notebooks/natcompsci2022/results_new/` | `pdb.id`, `pos.to` (0-based peptide index) | **derived/computed** — contact calls computed from **experimental** PDB crystal coordinates |
| `markup_2026.csv` | `notebooks/natcompsci2022/results_new/` | peptide length | **experimental** — sequences as deposited |
| `orient_metadata.json` | `src/tcren/data/` | `mhc.class` | **derived/computed** — tcren's per-structure annotation |

8,062 TCR↔peptide residue contacts over 370 structures, 19 (class, length) strata. Note what is
*not* used: no register call and no anchor call, class-I or class-II, so tcren's class-II register
heuristic — which is itself ported from mhcmatch — cannot make this circular.

Regenerate (from the `2026-mhcmatch-benchmark` repo):

    python bench/immuno/build_contact_profile.py \
        --tcren ~/vcs/code/tcren --out src/mhcmatch/data/contact_profile.py

**Validation gate (independent line).** The class-I 9-mer profile (176 structures, 3,772 contacts)
rank-correlates with Calis et al. 2013 (PMID 24204222) Table 2's KL position importances at
Spearman ρ = **0.943** (P3–P8, n = 6, one adjacent swap at P7/P8). Those KL values are derived from
immunogenicity *labels*; this profile is crystal *geometry* — the two share no data. P1/P2/P9 are
excluded because the paper itself marks them NA (anchors), which is also the conservative choice:
they are the three positions that would most inflate ρ. Caveat: only P4–P7 carry a significance star
in the source, so ρ = 0.943 leans partly on ranks the authors did not call significant.

## `complement_mhc1_human.json` / `complement_mhc1_mouse.json`

The parameters `mhcmatch.complement` scores with: a 30-feature design, its standardizer, fourteen
fitted log-odds tables, and a linear head. One file per species, **never pooled** — the mouse arm is
fitted on mouse rows only, because a pooled fit would let the larger human arm set the anchor and
TCR-face tables for both.

| part | what it is | provenance |
|---|---|---|
| `features` | the 30 column names, in the order `design()` builds them. **This list, not `complement.BLOCKS`, is the contract** — `_load()` fails on any length mismatch | **derived/computed** |
| `standardizer` | column means and standard deviations of the fitting arm | **derived/computed** |
| `log_odds` | thirteen 20-cell tables (amino-acid counts, by role and by length bin) and one 400-cell table (adjacent TCR-facing residue pairs) | **derived/computed** — fitted to experimental T-cell-assay labels |
| `logistic` | intercept and 30 coefficients on the standardized scale, ridge `tau = 4.0`. **This is what `score()` uses** | **derived/computed** |
| `fits.em`, `fits.supervised` | two class-conditional Gaussian fits, vendored so the comparison behind choosing the linear head stays re-checkable. **Not used at runtime** | **derived/computed** |
| `paratope` | TCRen contact potential marginalised over 28,250,990 TRB CDR3 loops | **derived/computed**, label-free |
| `anchors`, `alphabet`, `kd_threshold` | `(0, 1, 2, -2, -1)`, AA20, and the median Kyte–Doolittle value | **fixed** — `kd_threshold` is derived from the vendored `aa_tables` |

Fitted arms, one per file: `chowell_rebuilt/human` (464,161 rows, 14,712 immunogenic, prevalence
0.0317) and `chowell_rebuilt/mouse` (47,140 rows, 5,154 immunogenic, prevalence 0.1093), both from
`immunogenicity/chowell_rebuilt.tsv.gz` on the `isalgo/pmhc_data` dataset, seed 20260817. Positives
are peptides with a positive T-cell assay; negatives are eluted **self** ligands plus, for the human
arm, the HLA Ligand Atlas thymus immunopeptidome. The `aa` and `kmer` tables are refitted inside
every cross-validation fold so the reported AUROCs do not read their own labels back; the deposited
tables here are the whole-arm fit.

Regenerated in the benchmark repository, not here:

    python bench/neoag/complement.py --fit chowell_rebuilt --tables all

which writes `bench/neoag/complement_fit.json` and `complement_fit_mouse.json`; those are vendored
to these two paths on a version bump. Evidence — the block ablation, the four estimators, corpus
transfer and the size-matched cross-species comparison — is in `bench/results/complementarity.md`,
and the corpus construction rules and arm counts are in `bench/results/corpus_arms.md`. Which
features the model actually uses, including the nine Kidera factors it does **not** fit and an ESM2
comparison, is in `bench/results/complement_audit.md`. Nothing in this package recomputes any of it.

**Read the AUROCs as increments.** These corpora carry composition artefacts — cysteine marks the
label in the Chowell construction and reverses sign in the Kešmir one — so a 20-way composition
logistic alone reaches 0.68–0.74 on the Chowell arms under the same folds. `bench/results/
corpus_composition.md` and `KEY_FINDINGS.md` record the measurement.

## `complement.PARATOPE` (inline, not a data file)

The peptide-side reduction of the **TCRen** contact potential: `(mean, sd)` per peptide residue,
the receptor axis marginalised over a background repertoire. Read by `complement`'s `pot` block as
`para_tcr` / `para_sd_tcr`.

| input | where | what is taken | provenance |
|---|---|---|---|
| `TCRen_potential.csv` | `tcren` 2.8.0, `src/tcren/data/` | the **directed** 19×20 `J(CDR3 residue, peptide residue)` contact energy, never symmetrised and never reindexed to a square | **derived/computed** — a statistical contact potential fitted by tcren to solved TCR:pMHC interfaces |
| `human.trb.ntvj.vdjtools.tsv.gz` | `~/hf/airr_control` | the amino-acid composition of **28,250,990** TRB IMGT CDR3 loops, as the marginalisation measure | **experimental** — sequenced background repertoires |

- **Generated by** `bench/neoag/paratope.py` in
  [`2026-mhcmatch-code`](https://github.com/repseq/2026-mhcmatch-code) (private; released to reviewers);
  recorded in `bench/results/paratope_basis.md`.
- **Junction vs CDR3 is load-bearing.** `airr_control`'s `cdr3aa` column is the *junction* — it
  carries Cys104 and Phe/Trp118, so it is two residues longer than IMGT CDR3. Composing the measure
  from it unstripped would put a spurious cysteine in every receptor, and cysteine is precisely the
  residue TCRen leaves undefined on the receptor side.
- **Cysteine has no TCRen row**, so its repertoire share is redistributed over the residues the
  potential defines rather than being scored as zero energy.
- **The measure is the whole-loop composition, and that is now known to be the wrong one for a
  contact potential** — see `complement.PARATOPE_CONTACT` below, which fixes it. `PARATOPE` remains
  the default so every recorded number reproduces.

`MJ1996_contact_energies.csv` from the same `tcren` release backs the groove side of the same block
(`mj_anchor` / `mj_tcr`); the AAindex `MIYS850101` partition energy vendored into
[`aa_tables.py`](#aa_tablespy) is a different Miyazawa–Jernigan table and the two are not
interchangeable.

## `complement.PARATOPE_CONTACT` (inline, not a data file)

The peptide-side TCRen vector marginalised over the receptor residues that actually **contact**
peptide, offered beside `complement.PARATOPE` as `score(..., paratope="contact")`.

- **Generated by** `bench/immuno/paratope_contact_basis.py` in
  [`2026-mhcmatch-code`](https://github.com/repseq/2026-mhcmatch-code) (private; released to reviewers); result table
  `bench/results/paratope_contact_basis.md`, composition and cross-check in `paratope_contact.md`.
- **Two factors, each estimated where it can be.**
  `f_contact(a) = Σ_{clonotypes, i} P(contact | locus, i, L) · 1[cdr3_i = a]`. The geometry
  `P(contact | locus, i, L)` comes from tcren's 370 solved TCR:pMHC complexes (8,062 interface
  contacts); the residue identity from 28,250,990 TRB background clonotypes. Composition is
  **never** taken from the crystals — 370 complexes are a poor estimate of which residues sit
  mid-CDR3 and would inherit that set's TCR bias.
- **Why the shipped `PARATOPE` is the wrong measure for this potential.** TCRen is a *directed
  contact* potential, and only 35.4% of TRB loop residues ever contact a peptide. The germline
  V-encoded head and J-encoded tail that carry most of the flat whole-loop mass are exactly the
  positions the structures put at `P(contact) = 0.00` (TRB, length 12: positions 1–2 and 9–12), so
  conditioning on contact **is** the germline-flank trim and needs no separate rule.
- **Not a rescaling.** Spearman against `PARATOPE` is **+0.7549**, 19 of 20 residues change rank,
  and the spread widens from 0.2440 to 0.3255. A structure-free second route restricted to the
  N-D-N insert (located per clonotype by `VEnd`/`JStart`, no crystals at all) agrees at Spearman
  +0.7173, so the weighting is not an artifact of the structure set.
- **Pipeline check.** The same code path under the flat composition reproduces the shipped
  `PARATOPE` to max |diff| **1.10e-04**, so the two vectors differ by the measure and nothing else.
- **Cysteine has no TCRen row on the receptor side**, so its repertoire share is redistributed over
  the residues the potential defines rather than scored as zero energy. This vector is therefore
  not a second route to the cysteine artifact.
- **Opt-in.** The default stays `paratope="loop"`, so every recorded number for the shipped
  complementarity model reproduces unchanged.

## `corpus_tables.npz`

The **sliding-3-mer count tables** over each mimicry reference component's TCR face — the reference
side of the Łuksza-form corpus term that `C_corpus_thymus` / `C_corpus_self` / `C_corpus_viral`
read. One flat `20**3 = 8,000` float64 array per (class, component, species), plus a JSON `meta`
array carrying the build version and `k`.

**Derived, not experimental.** Every table is a deterministic function of deposits that ship in
`isalgo/pmhc_data` and the reference proteomes — no fitting, no parameters, nothing measured here.
Which deposit stands behind a key is `mimics.ref_path(category, species)`, and the three categories
resolve differently per species:

| component | human reference | mouse reference | how the species is selected |
|---|---|---|---|
| `thymus` | `thymus/thymus_immunopeptidome.tsv.gz` | `thymus/thymus_immunopeptidome_mmu.tsv.gz` | two files, `mimics.SPECIES_REFS` |
| `viral` | `ligandome/viral_foreign_iedb.tsv.gz` | the same file | one file, filtered on `mhc_species` |
| `self` | `proteome/human.fasta.gz` | `proteome/mouse.fasta.gz` | two files, `PROTEOME_REFS` |

| key | N (reference windows) | 3-mer cells filled, of 8,000 | build |
|---|--:|--:|--:|
| `mhc1\|thymus\|human\|3` | 140,482 | 6,083 (76.0 %) | 0.7 s |
| `mhc1\|thymus\|mouse\|3` | 25,264 | 3,877 (48.5 %) | — |
| `mhc1\|self\|human\|3` | 121,968,158 | 8,000 (100 %) | 51.4 s |
| `mhc1\|self\|mouse\|3` | 112,565,681 | 8,000 (100 %) | 35.3 s |
| `mhc1\|viral\|human\|3` | 136,618 | 7,340 (91.8 %) | 0.7 s |
| `mhc1\|viral\|mouse\|3` | 40,244 | 5,424 (67.8 %) | — |
| `mhc2\|thymus\|human\|3` | 1,996,006 | 7,311 (91.4 %) | 1.6 s |
| `mhc2\|thymus\|mouse\|3` | 186,758 | 5,011 (62.6 %) | — |
| `mhc2\|self\|human\|3` | 110,932,623 | 8,000 (100 %) | 14.5 s |
| `mhc2\|self\|mouse\|3` | 101,989,053 | 8,000 (100 %) | 10.0 s |
| `mhc2\|viral\|human\|3` | 1,205,107 | 7,866 (98.3 %) | 1.4 s |
| `mhc2\|viral\|mouse\|3` | 162,916 | 6,058 (75.7 %) | — |

The four mouse `thymus` and `viral` keys shipped in the artifact and appeared in no row of this
table until 2026-09-03; the build column is blank for them because its timings come from the build
in which only the eight human-and-self keys existed, and inventing a number for the others would be
worse than leaving it out. Rerun `mhcmatch build corpus` to fill them.

**Read the filled-cell column before reading a mouse corpus coefficient.** `self` is a proteome in
both species and fills every cell, so `C_corpus_self` is on the same footing in either host. The
two assayed channels are not: the mouse thymic deposit is **8,151 rows / 6,791 distinct peptides**
against the human file's 53,878, drawn from PRIDE PXD008733 (`H-2Db`, `H-2Kb`) and MassIVE
MSV000087031 (`I-Ab`), and the mouse slice of the viral file is **3,112 MHC-I rows against
51,972**. More than half of `mhc1|thymus|mouse|3` is therefore exactly zero, so a mouse
candidate's thymic density is read against a table with 4,123 empty cells and two H-2 allotypes
behind it. That is a property of what has been deposited in mouse, not of the construction, and it
is why the mouse-fitted corpus coefficients are reported nowhere.

**Why it is vendored: 145.4 kB of artifact against 115.6 s of per-process rebuild.** Regenerate on
a version bump, a deposit change or a change to what a face is:

```zsh
mhcmatch build corpus
```

The builder prints `** MOVED **` beside any table whose contents changed against the committed one,
so a rebuild that moves a number is visible at build time rather than in a downstream AUROC.
`tests/test_mimicry.py::test_the_vendored_corpus_tables_are_current_and_rebuild_bit_identically`
checks the version stamp, that every combination the builder declares is present, and bit-identity
against a live rebuild of the four deposit channels.

## `aggregate_mhc1.json` — the shipped neoantigen scorer (`EPIC`)

**Derived, not experimental.** A fitted model: nine standardised slopes in four hierarchical
blocks, one unpenalised intercept per screen and no global intercept, ridge `tau = 0.25`. The file
carries its own standardiser (`mu`, `sigma`) alongside the coefficients, because coefficients on
z-scores applied to raw axes move the *ranking* and not merely the calibration — the two travel as
one unit or not at all.

Read by `rank.aggregate()` and applied by `rank.aggregate_score()`. `mhcmatch rank --coefficients`
prints the terms and `--holdout` the per-screen held-out AUROCs; those two commands, not any
document, are the record of what a given install actually scores with.

**Two version vocabularies, told apart by shape.** `version` here is a *model* version and is an
**integer** — 3 through 11 — where a package version is dotted. `_build._stamp` therefore returns
`None` for this file and `mhcmatch build --check` **presence-checks it only**. It cannot tell a
current artifact from a stale one, which is deliberate (a model version does not move at every
release) and is why the copy below has to be checked by hand.

**Fitted in the benchmark repo, copied in by hand.** There is no in-process builder and no install
script; `_build.EXTERNAL["aggregate"]` prints the command and stops.

    # in ~/vcs/projects/2026-mhcmatch-benchmark, and only after `corpus` and `features`:
    python bench/epic/fit.py --physchem rose_af5 --presentation binder --density log10a
    # then, deliberately:
    cp bench/epic/aggregate_mhc1.json ~/vcs/code/mhcmatch/src/mhcmatch/data/aggregate_mhc1.json

`fit.py` writes the candidate and logs that it was **not** copied. Diff the two files before and
after the copy: nothing else will catch a bad one.

**Corpus.** `bench/epic/neoantigens.parquet` + `bench/epic/features.parquet`, both built from the
`isalgo/pmhc_data` deposit and both stamped with the `mhcmatch` version that wrote them —
`bench/epic/optimize.py` refuses any other, so a fit cannot silently read a frame built by a
different library. The screens the fit spans, its row and positive counts and its BIC are all
recorded in the artifact's own `fit` block.

**History.** v3 fitted `binder`; v4 respecified it to `pres` on the argument that `occupancy` is a
monotone function of the same predicted Kd; v5 refitted v4's specification on a rebuilt corpus; v6
returns to `binder`, because that argument conflated a within-allele `%rank` with an absolute Kd —
measured, `pres` is *more* collinear with `binder` (Spearman +0.8797) than `occupancy` is (+0.7431).
`bench/results/epic_binder_vs_pres.md`. v7 puts the density term on its natural scale, `log10a`
in place of `occupancy`; v8 fits abundance as a level rather than a rank; v9 splits expression into
the two terms `expr_lvl` and `expr_norm`, reaching the nine that ship.

**v10 (2026-08-28) is v9's nine terms refitted, and it is the first entry here that ships against
its own verdict.** Nothing about the specification moved — same features, same blocks, same
`tau = 0.25`. What moved is underneath: `mhcmatch.affinity.PottsAffinity.predict_y` gained the
corpus length prior and `mhcmatch.calibrate.percent_rank` gained an extrapolated upper tail, so
`binder` and `log10a` are computed from a binding layer that knows two things v9's did not. Only
those two standardisers move; the other seven `mu`/`sigma` pairs are bit-identical to v9's, which
is the check that the refit is downstream of the binding change and of nothing else.

    BIC        4390.2 -> 4328.3
    LOO mean   0.6942 -> 0.6998   (leave-one-screen-out, 8 screens)
    coef       log10a +0.2914 -> +0.4005 (+37%), binder +0.5481 -> +0.4623 (-16%)

All nine terms stay individually significant (|z| 2.52-6.22, max p = 0.0118, sign stability >= 0.995
over 400 clusters). The fit leans harder on affinity because affinity now knows something.

**v11 (2026-08-29) supersedes a v10 that no longer reproduces, and is fitted on a corpus that
changed twice.** The specification is untouched again — same nine features, same blocks, same
`tau = 0.25`. Three things moved, and they are separable:

*1. A v10 that does not rebuild.* Every fit the chain produces gives `binder` `(mu, sigma)`
`(-1.392931, 0.514058)` and `log10a` `(-2.865114, 0.924722)`; the shipped v10 file carries
`(-1.437181, 0.505559)` and `(-2.878300, 0.915829)` — on exactly the two columns the v10 note above
says define it over v9, and on no others. The intent recorded there is still the right description
of what v10 was *for*; the file drifted from it. The cause was not established. On the author's
decision the artifact is **superseded rather than reconciled**, which restores the property that a
shipped artifact rebuilds from the chain.

*2. Parent genes resolved.* `expr_lvl` and `expr_norm` key on a gene symbol, and the symbol was
missing on 356,387 of 695,811 corpus rows (51.2%) and on 5,205 of 5,833 positives (89.2%). Every
such row took one mean-imputed value, so on VACCIMEL `expr_norm` had standard deviation **exactly
0.0000** and AUROC **exactly 0.5000** while carrying v10's second-largest coefficient. Resolved by
`bench/neoag/resolve_genes.py` (seqtree, radius 2 — a neoantigen can carry more than one mutation);
coverage 48.8% -> 99.5%. `bench/results/gene_resolution.md`.

*3. Gfeller_GBM removed from the fit.* It is not an independent screen but the compendium published
with Koohy's immunogenicity comparison, and **2,733 of its 2,833 (peptide, allele) pairs (96.5%)
are Gfeller**, which `corpus.py` rule 10 already excludes as viral and self rather than neoantigen.
It also carried 592 CEDAR pairs, 100 of GBM's 150 and 49 of ITSNdb's 197, so leave-one-screen-out on
Gfeller left 83.4% of Gfeller's own pairs in training and on GBM 66.7%. Removing it costs 144 of 741
positives — 19.4% — and is recorded as `corpus.py` rule 10a.

    rows       342,432 -> 339,599      positives  741 -> 597      screens  8 -> 7
    LOO mean   0.6998 -> 0.7102        (leave-one-screen-out)
    BIC        4328.3 -> 3109.8        NOT comparable: n and the positive count both moved
    coef       binder +0.4623 -> +0.7569 (+64%), expr_norm +0.4950 -> +0.2155 (-56%),
               expr_lvl +0.3704 -> +0.5180, log10a +0.4005 -> +0.1713

**The two large coefficient moves have one mechanism between them.** 2,733 viral-and-self pairs were
pulling the fit away from presentation and onto the corpus and expression channels; removing them
lets `binder` carry what it actually carries, and giving `expr_norm` real variation drops it to what
the information supports. Held-out accuracy is unchanged within noise despite losing a fifth of the
positives, which says those rows were bending the fit rather than carrying it.

**One number to read with care:** the bootstrap resamples (patient, screen) clusters, and the count
fell 3,294 -> 527 with Gfeller_GBM, whose rows each formed their own cluster. Intervals are wider
than v10's for that reason and not because the fit is less certain.


**Its `verdict` block reads `"ship": false` and it ships anyway — that is the author's call, taken
2026-08-28, and it is recorded rather than edited away.** `fit.py`'s bar is *no per-screen
regression at all*, and v10 posts 3 improvements, 3 ties and 2 regressions: IEDB_neoag
0.7270 -> 0.7217 (-0.0053, 235 discordant pairs, against that screen's own 1/234 = 0.0043
resolution) and ITSNdb 0.5510 -> 0.5329 (-0.0180, 159 discordant pairs of 8,832). Both are on
screens where the measurement is thin -- ITSNdb sits near chance under *both* artifacts, and it
also blocks the currently-shipped v9 by the same bar. Against them stand NCI +0.0134 over 449,998
discordant pairs, TESLA +0.0320 and Gfeller_GBM +0.0285. The bar is not relaxed and the verdict is
not rewritten; the artifact carries its own dissent.

## `aggregate_mhc1_mouse.json` — the mouse class-I scorer

**Derived, not experimental**, and the same object as the human artifact in every structural
respect: nine standardised slopes under the same nine feature names in the same order, ridge
`tau = 0.25`, its own `mu`/`sigma`, a *model* version that is an integer (**`3`**) and a `release`
that is dotted (`1.12.0` -- the package version the fit was **accepted** in, which is what a
manuscript cites). `rank.aggregate(cls, species, mode)` resolves it through
`rank.AGGREGATE_ARTIFACTS`, keyed `(cls, species, mode)`. All four `(cls, species)` cells are
fitted from 1.12.0; `pathogen` remains a registered mode with no shipped artifact rather than a
silent alias for `neoantigen`, and an unregistered species raises instead of being served a
neighbour's coefficients.

The class-II half of this section moved to **`aggregate_mhc2_human.json` /
`aggregate_mhc2_mouse.json`** below: from 1.12.0 the two class-II fits are one specification with
six terms and no corpus block, and describing them beside a nine-term class-I fit made the shared
feature list look like an accident.

**Nine coefficients, seven free parameters, since model version 2.** (Class II is six and six --
see below.) The corpus block is one fitted scalar on the *human* artifact's corpus direction (`+0.3255, -0.8601, +0.3928` after normalising
v11's `C_corpus_*`), not three independent coefficients -- a projection with one free scalar **is**
three coefficients constrained to be proportional. The file lists all nine so one library scores
both species with no branch. Every per-term array (`coef`, `sd`, `boot_sd`, `z`, `p`, `ci95`,
`sign_stability`, `mu`, `sigma`) is expanded to nine alongside them;
`test_every_registered_artifact_declares_the_features_it_carries_coefficients_for` is the guard,
added after the first v2 copy shipped five of them at seven and made
`rank --coefficients --species mouse` raise `IndexError` on the eighth row.

    # in ~/vcs/projects/2026-mhcmatch-benchmark-ext:
    MHCMATCH_MODEL_RELEASE=1.12.0 python bench/pmhc_data/clean_neoantigens.py
    MHCMATCH_MODEL_RELEASE=1.12.0 python bench/epic/fit_mouse.py --cls mhc1 --model-version 3
    # `--corpus-axis human` is the default and is what version 2 means; `--corpus-axis free`
    # reproduces version 1's three independent corpus coefficients.
    # then, deliberately:
    cp bench/epic/aggregate_mhc1_mouse.json ~/vcs/code/mhcmatch/src/mhcmatch/data/

**Corpus.** `~/hf/pmhc_data/neoantigens/neoag_tested_mmu.tsv.gz`, the IEDB mouse neoantigen
deposit, at class I's own peptide lengths and keyed on `mhc_a_pred`: **921 rows, 379 immunogenic,
61 references, 6 H-2 allotypes.** (**923 / 380 through version 2**: the deposit was cleaned of
pathogen and unattributable rows for 1.12.0 and two of them were in this fit. Nine terms and the
pinned corpus axis are unchanged; the version moves because a citation has to name one fit.) No row is dropped for a missing term: seven of the nine are
populated on every row and the two expression terms are imputed to the population median, which is
the convention `rank.aggregate_score` documents.

Every column is computed by `mhcmatch rank pairs --species mouse --score features`, so the panel,
the three corpus channels and both expression terms are the library's own **mouse** references
throughout. That flag is load-bearing rather than cosmetic: the human panel holds no `H-2-Kb`
ligands, so the groove is borrowed from its human kernel neighbours and SIINFEKL -- the canonical
H-2-Kb binder -- scores `binder` **-1.148** against **+1.523** under `--species mouse`.

**One unpenalised intercept per reference, and the fit is not interpretable without it.** The human
artifact gives each *screen* an intercept so prevalence and candidate generation stay out of the
slopes. This deposit is one screen and 61 publications whose positive rate runs 0 % to 90 %, so the
reference is where that variation lives. Against a single pooled intercept every one of the nine
coefficients came out at or below zero and the held-out figure sat at 0.4633.

**No held-out split is fitted or shipped, on any of the three artifacts.** These are GLMs: the
deliverable is a coefficient and the interval around it, and the interval is a **cluster bootstrap
over `reference_id`** -- the publication is the unit that repeats in these deposits, so that is
what is resampled. No artifact carries a `cv_*` block, and `bench/epic/fit_mouse.py --folds 0`
(the setting all three shipped fits were run at) fits none. Where a discrimination figure does
appear it is **in-sample AUROC within a reference**, macro-averaged over the references carrying
at least three of each class, and it is a fit diagnostic rather than a claim: a pooled figure over
this deposit compares a paper reporting 144 positives of 201 against one reporting 1 of 191, so it
would be mostly a base-rate difference between laboratories.

**One thing to know before quoting a coefficient.**

*`expr_norm` is negative in mouse (-0.2314, z -2.78) where it is positive in human (+0.2155), and
the two are not measuring the same thing.* There is no mouse tumour transcriptome: a tumour
abundance exists for 34 % of rows (GEO GSE281579) and for the rest `rank._expression_for` falls
back to the gene's FANTOM5 tissue median, which is what `expr_norm` already is -- so on those rows
the two columns are identical. The sign was checked against that and survives it: under the `pan`
arm, where `expr_norm` is the gene's pan-tissue median and by construction never equals `expr_lvl`,
the coefficient is **-0.2432 at z -3.31**, larger and more significant. The indicator itself
(`expr_observed`, +0.3136, z +0.75) is null, so which publication deposited an abundance is not
carrying it either. `bench/results/epic_mouse_fit_mhc{1,2}.md` has all four arms.

## `aggregate_mhc2_human.json` / `aggregate_mhc2_mouse.json` — the class-II scorers

**Derived.** One specification, two species: **six standardised slopes** under the same six names
in the same order -- `binder`, `log10a`, `expr_lvl`, `expr_norm`, `C_phys_buried`, `C_phys_charge`
-- ridge `tau = 0.25`, each with its own `mu`/`sigma`, one unpenalised intercept per reference, and
no intercept in the shipped file. `rank.TERMS_MHC2_EXPECTED` names the six;
`test_both_class_II_artifacts_carry_the_same_six_terms_and_no_corpus_block` is the guard that keeps
the two comparable term by term. Human is model version **`1`**, mouse **`3`**, both `release`
`1.12.0`.

**There is no corpus block, and that is a statement about the reference and not about the
coefficient.** A `C_corpus_*` channel is a Łuksza density over a reference set of peptides --
thymic, self, viral -- and what is deposited for all three is a **class-I** set. Contracting a
15-mer class-II register against a 9-mer density is not a weak feature, it is the wrong question,
so the block leaves the design entirely: no `C_corpus_*` column is fitted, `blocks` lists three
entries rather than four, and the four corpus-geometry keys (`corpus_k`, `corpus_mask`,
`corpus_kernel`, `corpus_shapes`) are **absent** rather than declared-and-unused. The nine names
stay in `rank.AGGREGATE_FEATURES`, because a recorded result cites the model version that produced
it and a registry that drops a name cannot say what those numbers were.

    # in ~/vcs/projects/2026-mhcmatch-benchmark-ext:
    MHCMATCH_MODEL_RELEASE=1.12.0 python bench/pmhc_data/clean_neoantigens.py
    MHCMATCH_MODEL_RELEASE=1.12.0 python bench/mhc2_human/fit_human_mhc2.py
    MHCMATCH_MODEL_RELEASE=1.12.0 python bench/epic/fit_mouse.py --cls mhc2 \
        --corpus-axis none --model-version 3
    # then, deliberately:
    cp bench/epic/aggregate_mhc2_human.json           ~/vcs/code/mhcmatch/src/mhcmatch/data/
    cp bench/epic/aggregate_mhc2_mouse_axis-none.json \
       ~/vcs/code/mhcmatch/src/mhcmatch/data/aggregate_mhc2_mouse.json

**Corpus, human.** `~/hf/pmhc_data/neoantigens/cedar_neoag_mhc2_hsa.tsv.gz` -- **1,112 rows, 656
immunogenic (59.0 %), 157 references, 72 allotypes, 238 genes.** Built from the CEDAR export by
`bench/mhc2_human/build_cedar_human_mhc2.py`, one row per `(peptide, allele)`, label positive if
any assay on that pair read `Positive*` and negative only when every one read `Negative` (149 pairs
disagreed and took the positive). Every antigen is a human self protein: the pathogen rows are a
different mode and are filtered out, which is also what took gene coverage from 43.6 % to **98.2 %**
-- a viral peptide has no host gene, so the missingness was never a resolution problem.

**This deposit is human self-antigen CD4 response, and it is not a tumour cohort.** 364 of the
1,112 rows are type 1 diabetes, 260 are healthy donors, 63 rheumatoid arthritis, 51 multiple
sclerosis; every cancer together is **143 rows**. The mechanism -- a CD4 response to a self peptide
-- is the one the model is about, and a consumer scoring tumour neoantigens with it is
extrapolating. The composition is in the artifact under `fit.population` so that can be read
without the report, and `bench/results/epic_human_fit_mhc2.md` has all 46 diseases.

**`expr_norm` needs a tissue or it is a copy of `expr_lvl`.** With no `--tissue`, `rank` gives both
terms the gene's pan-tissue median over the same floor and the design carries one column twice. The
disease resolves the context, library first: `expression.resolve_context` already maps 21 of
CEDAR's disease strings -- every cancer it names -- onto TCGA codes and their matched GTEx normals,
and `fit_human_mhc2.DISEASE_TISSUE` adds 11 autoimmune target organs it cannot know (T1D to
pancreas, MS to cortex, pemphigus to skin). **686 of 1,112 rows (61.7 %) over 30 contexts**; the
rest take the pooled human floor and `expr_norm` falls back to the gene's pan-tissue median. RA,
lupus and Vogt-Koyanagi-Harada are deliberately unmapped: GTEx samples no synovium and no uvea.

**What the human class-II fit says.** BIC **1595.8** on 1,112 rows, deviance 452.5, 157 intercepts
at 157.0 effective df. Two of the six coefficients are sign-stable across the 400-resample cluster
bootstrap on `reference_id` -- `binder` **+0.3773** at 0.98 and `C_phys_buried` **+0.1710** at 0.97
-- and the other four are not, `log10a` at 0.63, `expr_lvl` 0.65, `expr_norm` 0.61,
`C_phys_charge` 0.81. Alone and in sample, within reference: `binder` **0.5645**, `C_phys_buried`
**0.5485**, `log10a` **0.5386**, and **both expression terms below chance** (`expr_lvl` 0.4035,
`expr_norm` 0.4393) -- which is what an autoimmune-dominated deposit should look like, since an
islet antigen's abundance is not a neoantigen's. In sample within reference the joint fit reaches
0.6020.

**Corpus, mouse.** Unchanged from version 2 -- `neoag_tested_mmu.tsv.gz` unioned with CEDAR's
mouse-derived H-2 assays, **468 rows, 177 immunogenic, 30 references, 7 allotypes**. Only the
corpus block moved. Arm-vs-arm on that same population, `vanilla`, v2 → v3: BIC **562.3 → 556.2**,
which is `log 468 = 6.15` -- one parameter's worth, and the block was spending three names on it.
`bench/results/epic_mouse_fit_mhc2_axis-none.md`.


## `anchor_model_*_mouse_*.pkl.gz` — the mouse AnchorModel pickles

**Derived.** Five files, the mouse counterpart of the human registry entry for entry, built by
`mhcmatch build anchor` from `Store.from_pmhc(species="mouse")` and loaded read-only under the same
version / `panel_sha` / params guard. 0.66 MB between them against the human set's 5.3 MB, and 12 s
to build against 304 s, because the mouse panel is a twentieth of the human one.

They are a performance artifact and nothing else. `panel_sha` already made a mouse run *correct* --
a mouse panel misses every human entry and falls through to a fit -- so what these change is that
it no longer refits on every call, including the class-II register+EM the pickles exist to avoid.
`diffusion.vendored_models(species)` selects the registry; `Store.species`, set by `from_pmhc`, is
how the loader knows which to try.

## `known_neoantigens.idx` / `known_neoantigens.json` — the validated-epitope whitelist index

**What it is.** A `seqtree.Index` over every peptide `mhcmatch.known.load()` collects into its
`neoantigen` set: 23,299 sequences, lengths 4 to 50, each one an assay called **positive**. It backs
`--keep-epitopes builtin`, the whitelist of candidates with a validated T-cell response that no
`--rank-threshold` may drop. The `.json` is a version sidecar; a `seqtree` index is opaque binary and
cannot carry a stamp, so `build --check` reads the sidecar beside it.

**Source.** Five deposits, listed in `known.SOURCES["neoantigen"]` and repeated in the sidecar, all
fetched from `isalgo/pmhc_data`:

| deposit | what it contributes |
|---|---|
| `neoantigens/neoag_tested.tsv.gz` | the aggregated cohorts (CEDAR, Gfeller, Neopep, ITSNdb, TESLA, GBM, VACCIMEL), positive rows |
| `neoantigens/neoantigens_tested_peptides.tsv.gz` | the epitope-resolution screens (NCI, HiTIDE, TESLA), positive rows |
| `neoantigens/nci_gartner_mmp.tsv.gz` | the NCI/Gartner deconvolved minimal peptides |
| `neoantigens/neoag_tested_hsa.tsv.gz` | IEDB neoantigen records, human |
| `neoantigens/neoag_tested_mmu.tsv.gz` | IEDB neoantigen records, mouse |

**Regenerate.** `mhcmatch build known` — 1.3 s, of which 1.2 s is the deposit scan and 0.02 s the
index build and save.

**Why it ships rather than being assembled on demand.** The five deposits carry roughly 950,000 rows
between them, so building the set is a download plus a full-file scan. A thousand-sample Nextflow run
would pay that a thousand times, or race on whatever cache it wrote to avoid doing so. Pre-built, the
index reloads in **~1 ms** and answers **~1.45 M queries/s** through `search_batch`, which releases
the GIL and uses every core. Nothing on the predict path ever builds it.

**Provenance is derived, not experimental**: the peptides are experimental (each is an assay result),
the index over them is computed. It adds no peptide that `known.load` does not already return, so it
is an index and never a second definition of the corpus.
