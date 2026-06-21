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
