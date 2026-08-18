# `mhcmatch.vector` at v0.14.0 — what it does, and where the evidence says it is thin

Audit date 2026-08-18, against `src/mhcmatch/vector.py` (934 lines) and `tests/test_vector.py`
(683 lines). Evidence tiers and citations are in [`vector_evidence.md`](vector_evidence.md).

## What is shipped

| # | capability | entry point | basis |
|---|---|---|---|
| 1 | withdraw unsafe units before spending capacity on them | `screen`, `self_origin_risk` | MAGE-A3/titin and MAGE-A12/brain fatalities |
| 2 | choose how many units, per allotype | `select` | per-allotype saturating objective |
| 3 | lay them out and pick a spacer | `order`, `scan_junctions`, `junction_windows` | pVACvector's junction objective |
| 4 | build a unit from context | `unit`, `units_from_context` | 27-mer centred, Kreiter/Bijker |
| 5 | emit and repair nucleotides | `back_translate`, `slippery_sites`, `deslip`, `translate` | Mulroney 2023 |
| 6 | plumbing | `store_binder`, `rebuild`, `from_sequence`, `Unit`/`Selection`/`Cassette` | — |
| 7 | one-call pipeline | `mhcmatch vector`, `mhcmatch deslip` | — |

Three design properties worth keeping, because they are what make the module auditable:

- **`order` and `scan_junctions` take an injected `binder` callable**, so layout is testable with no
  `Store`, no panel and no download.
- **`_greedy_2opt` uses no RNG**, so a cassette is a pure function of its inputs — unlike
  pVACvector's simulated annealing.
- **`select` refuses to default `n0`.** The literature does not fix it, so a caller who has not
  chosen has not finished designing.

## Gaps

Ordered by how much the retrieved evidence says they cost.

### G1 — A mixed CD8 + CD4 cassette cannot be laid out *(highest)*

The register vocabulary is **one tuple per cassette, chosen once**. `order()` and
`scan_junctions()` take a single `lengths` argument, and `cli.py:862` resolves it from a single
`--cls` flag — `JUNCTION_LENGTHS = (8,9,10,11)` for `mhc1`, `MHC2_JUNCTION_LENGTHS = (12,13,14,15)`
for `mhc2`. The CLI comment states the design directly: *"One register vocabulary for the whole
command."*

That is correct for a single-class cassette and cannot express a mixed one. Run it as `mhc1` and
every class-II junctional epitope — the exact hazard Livingston 2002 measured — is invisible; run it
as `mhc2` and the class-I registers are. `select(cls=...)` likewise filters to one class, so the two
payloads can only be selected separately.

This matters because linkage is not optional: co-delivering the CD4 and CD8 components as **separate
constructs produced no antitumour immunity at all**, where the fusion worked (PMID 15270727). The
module can currently build either half and not the thing the evidence says to build.

### G2 — Spacer choice is not class-conditional

`SPACERS` is one global tuple in one fixed order. The evidence splits cleanly by class: the only
head-to-head MHC-I processing experiment favours **alanine-based** spacers over `GGGS`
(PMID 36820900), while every `GPGPG` rescue result is class II or antibody (PMID 12023344,
PMID 22922658).

The module docstring currently argues for Gly/Pro-rich spacers from ligand-flanking *composition*
(PMID 30645615) — a bioinformatic argument — where a presentation assay went the other way for
class I. The ranking is defensible per class and misleading as one global list.

### G3 — Nothing models the N-terminal flank

`unit()` centres the mutation and stops. But TAP prefers N-terminally **extended** precursors, and
several real epitopes are poor TAP substrates as minimals and transported far better with up to four
extra N-terminal residues (PMID 9764810). Flanking effects are large and sometimes absolute — one
added C-terminal methionine restored presentation that was otherwise zero (PMID 9029109); two
flanking prolines maximised processing (PMID 26018465); a substitution five residues upstream
abolished an epitope without touching it (PMID 26446603).

A 27-mer gets native flanks by construction, which is why the default is sound. But the **spacer
replaces the native flank on one side of every junction**, and nothing asks whether the replacement
is good. PolyCTLDesigner (PMID 24107711) does exactly this and is the prior art we do not yet match.

### G4 — `order()` optimises only junctional binding

The edge cost is the strongest predicted binder spanning a junction. That is pVACvector's objective
and it is the right *first* one. It says nothing about whether the chosen spacer liberates its two
neighbours efficiently — which is the axis the one direct experiment actually measured
(PMID 36820900). A junction can be binder-clean and still suppress presentation of both units beside
it.

### G5 — There is no helper layer, and no per-unit help flag

CD4 help-dependence is **per-epitope**: of three epitopes in one antigen, only the immunodominant one
collapsed without help, and adding helper epitopes enhanced only that one (PMID 21810614). There is
no `Unit` field for it, no universal-helper slot (PADRE has both preclinical and clinical-stage
support: PMID 32145473, PMID 27903079), and no way to express "this unit needs help and that one
does not".

### G6 — No backbone

`back_translate` emits a bare CDS. Cap, 5'UTR, Kozak, signal peptide, MITD, stop, 3'UTR and poly(A)
are all outside the module, so the artifact it emits is not the thing that gets synthesised. Note the
caution before adding tags with a promised mechanism: a head-to-head of tPA, ubiquitin and LAMP-1
found all three beat the untagged vaccine but **none steered the arm it was chosen for**
(PMID 19356616).

### G7 — The layout solver explores one path through a large free set

Junction-free orderings are "astronomically" abundant (PMID 20033850), and no retrieved experiment
distinguishes them — position had minimal impact in the one assay that measured it (PMID 36820900),
and all six permutations of a three-region construct gave no junctional response (PMID 7521933).

`_greedy_2opt` returns one cheap path. That is the correct amount of effort for the *first* objective
and it wastes the freedom left over. The opportunity is a **second** objective inside the
junction-clean set, not a better search for the first.

## Not a gap

- **Duplication.** A centred 27-mer already contains every 8–14mer register spanning the mutation, so
  duplicating for register or allotype coverage buys nothing. Where duplication is genuinely forced,
  the requirement is a flexible separator — a tandem repeat was processed **only** when the copies
  were separated by two glycines (PMID 14993639). Worth supporting as a constrained option (G4 in the
  roadmap), not as a default.
- **Nesting.** Nested helper-over-CTL works, but distant help worked as well and position inside the
  nest did not matter (PMID 21810614). Linkage is what is required; geometry is not.
- **Frameshift repair.** `slippery_sites`/`deslip` already implement the Mulroney 2023 fix
  (PMID 38057663).
