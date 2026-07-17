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
(see `appendix/mhcmatch.tex` §4). Regenerate with:

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
