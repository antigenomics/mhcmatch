# Vendored data provenance

## `mhci_pseudo.fa` / `mhcii_pseudo.fa`

NetMHCpan-style **34-residue MHC pseudosequences** — the polymorphic groove positions that contact
the peptide (class I: α1/α2 of the MHC heavy chain; class II: α1 of the α-chain + β1 of the
β-chain). Header format `>ALLELE|n=<count>`; `X` marks an ambiguous/unresolved position.

- **MHC-I:** 4143 alleles (human HLA-A/B/C, mouse H-2, and other species).
- **MHC-II:** 2209 alleles (human HLA-DR/DQ/DP, mouse H-2 I-A/I-E, others).

Copied verbatim from the sibling `antigenomics/tcren` repository
(`tcren-ms/src/tcren/data/{mhci,mhcii}_pseudo.fa`, built by its `scripts/build_pseudo_fasta.py`),
which derives them from the NetMHCpan pseudosequence definition. Used by `mhcmatch.pseudoseq` as the
allele-similarity alphabet for the cross-allele diffusion model (see `appendix/mhcmatch.tex` §4).

These files are static reference data and small (~340 KB total), so they are vendored rather than
fetched. Re-sync from `tcren` if the pseudosequence definition is updated upstream.

## `structural_pockets_mhc1.tsv` / `structural_pockets_mhc2.tsv`

Per-anchor **structural pocket weights**: for each peptide anchor (MHC-I P1/P2/P3/PΩ-1/PΩ; MHC-II
P1/P4/P6/P9), the frequency with which each of the 34 groove pseudosequence positions makes a
heavy-atom contact (<5 Å) with that anchor residue, measured over pMHC crystal structures. Used as a
data-independent alternative to the learned-MI groove weights via `AnchorModel(weights="structural")`.

Measured by `bench/structural_pockets.py` from the **Canonical2026** TCR:pMHC structure set
(`antigenomics/tcren`, 372 usable structures): the 34-mer pseudosequence is threaded onto each groove
with tcren's C++ fitting aligner (`tcren._align`; no mmseqs/arda). Class per structure is assigned by
best pseudosequence fit (MHC-I single chain vs MHC-II α1+β1 chain-pair), giving **279 MHC-I** and
**93 MHC-II** structures. Regenerate with:

    conda run -n tcren-nb python bench/structural_pockets.py \
        --structures ../tcren-ms/data/Canonical2026 --out src/mhcmatch/data

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
