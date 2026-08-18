# `mhcmatch.vector` — development roadmap

Companion documents: [`vector_audit.md`](vector_audit.md) (what is shipped, where it is thin) and
[`vector_evidence.md`](vector_evidence.md) (the PubMed scan of 2026-08-18, every claim tiered).
Gap labels `G1`–`G7` below are the audit's.

## The organising claim

The four assembly questions have different amounts of evidence behind them, and the module should
spend its complexity accordingly:

| question | state of evidence | what that implies |
|---|---|---|
| what to **refuse** | two fatal precedents, mechanism known | hard exclusion — shipped |
| **how many** | no objective function published anywhere | one declared free parameter, `n0` — shipped |
| **what between** | one head-to-head MHC-I assay, one class-II rescue | class-conditional defaults — **V1** |
| **in what order** | three studies say position barely matters | satisfy the constraint, then stop — **V4** |

The single most consequential finding of the scan: **ordering is a constraint-satisfaction problem,
not an optimisation problem.** Junction-free layouts are astronomically abundant, and no retrieved
experiment distinguishes them. Effort spent on a better path search is effort not spent on the axes
that measurably move presentation — spacer identity and flanking sequence.

The second: **CD4 and CD8 payloads belong in one molecule.** Co-delivering them as separate
constructs produced no antitumour immunity where the fusion worked. That closes a fork the earlier
design memo left open, and it makes G1 the first thing to build.

---

## V1 — Class-aware assembly *(closes G1, G2)*

The one release that changes what the module can produce rather than how well it does it.

**V1.1 Per-junction register vocabulary.** Today the vocabulary is fixed once per cassette from a
single `--cls` flag (`cli.py:862`). Make the register set a function of the two units flanking each
junction instead: class-I ↔ class-I scans `(8,9,10,11)`, anything touching a class-II unit also scans
`MHC2_JUNCTION_LENGTHS`. Both tuples already exist; what is missing is resolving them per junction
rather than per command.

**V1.2 Per-unit binder alleles.** A class-II junction must be scored against the recipient's DR/DP/DQ
allotypes, not their A/B/C. `store_binder` takes `cls` at construction, so the natural shape is a
`{cls: binder}` mapping resolved per junction.

**V1.3 Class-conditional spacer defaults.** Replace the single `SPACERS` tuple with a selection keyed
on the junction's class pair:

- class-I ↔ class-I → alanine-based first (`None`, `AAA`, `AAY`), on the only head-to-head MHC-I
  processing assay (PMID 36820900)
- any junction touching class II → `GPGPG` first, on the only causal junctional-suppression rescue
  (PMID 12023344, replicated PMID 22922658)

`None` stays first in both: a clean junction needs no spacer, and every inserted residue is translated
sequence that could itself bind.

**V1.4 Mixed-class `select`.** `select(cls=None)` already groups by allele, which separates the
classes incidentally because the allele strings differ. Make that explicit and give class II its own
`n0` — class-II capacity is not class-I capacity and pretending otherwise hides an assumption.

**Verify.** Rebuild a cassette from a known payload with and without class-II units and confirm the
class-II junction registers are scanned and scored. Regression: a class-I-only payload must produce
the byte-identical cassette it produces today.

**Correction this release must carry.** The module docstring argues Gly/Pro-rich spacers sit in the
"permissive zone" from ligand-flanking composition (PMID 30645615). The one presentation assay went
the other way for class I. Restate it as class-conditional and name the assay.

---

## V2 — Flanking and processing *(closes G3, G4)*

**V2.1 TAP-aware N-terminal extension.** Score each unit's N-terminal residues against the human TAP
binding motif and, where a unit begins with residues deleterious to transport, extend it N-terminally
by up to four residues of native context — the extension length that raised ER transport measurably
(PMID 9764810). Only ever extend into **native** sequence, never invent it: an invented flank is a new
junction with none of the safety screening.

**V2.2 A liberation term in the junction objective.** `order`'s edge cost is currently the strongest
predicted binder spanning the junction. Add a second, explicitly-weighted term for whether the spacer
supports cleavage that liberates its two neighbours. Keep it opt-in and keep the weight a named
argument — the evidence establishes that spacers change liberation efficiency (PMID 36820900), not by
how much relative to junctional binding.

**V2.3 Report the flank each unit actually got.** `Unit` records `mutation_index`; it should also
record how much native context sits either side, because near a protein terminus the window clamps
and the flank silently shortens.

**Verify.** The TAP extension must be a no-op on units that already start with a favourable residue,
and every extension must be recoverable from the source context. Both objectives must be reportable
side by side on the same payload, because — as with `sum` vs `rate` — they will sometimes disagree,
and a caller who has not chosen has not finished designing.

---

## V3 — The helper layer *(closes G5)*

**V3.1 A help-dependence field on `Unit`.** Help-dependence is per-epitope, not per-cassette: of three
epitopes in one antigen only the immunodominant one collapsed without help (PMID 21810614). Add a
tri-state — needs help / does not / unknown — defaulting to unknown, and let it be set from evidence
rather than predicted, since nothing retrieved predicts it.

**V3.2 A universal-helper slot.** Allow one or more helper units that are not neoantigens — PADRE is
the obvious first, with preclinical CD4+CD8 gains (PMID 32145473, PMID 30999007) and clinical-stage
use as a genetic fusion without autoimmunity (PMID 27903079). It is furniture, so it should be
declared, not selected: it does not compete for `select`'s per-allotype budget.

**V3.3 Placement is free, linkage is not.** Do **not** implement nesting geometry. Distant help worked
as well as nested help and position inside the nest did not matter (PMID 21810614). What must be
enforced is that helper and helped units end up in the **same cassette** — which is the one thing
shown to matter, since separate constructs failed outright (PMID 15270727).

**V3.4 Constrained duplication.** Support duplicating a unit only with a mandatory flexible separator:
a tandem repeat was processed only when the copies were separated by two glycines (PMID 14993639).
Default off, and document that a centred 27-mer already carries every register spanning the mutation,
so duplication is for tiling frameshift/fusion ORFs, not for coverage.

---

## V4 — Layout freedom and the backbone *(closes G6, G7)*

**V4.1 Enumerate the junction-clean set, then choose within it.** Junction-free orderings are
astronomically abundant (PMID 20033850). Rather than returning one greedy path, return the clean set
(or a bounded sample of it) and apply a **second** criterion inside it — e.g. spread same-allotype
units apart, or put help-dependent units near the helper. State the second criterion as an assumption,
because no retrieved experiment ranks clean layouts.

**V4.2 Backbone assembly.** Give `back_translate` an optional backbone: cap, 5'UTR, Kozak, signal
peptide, MITD, stop, 3'UTR, poly(A). Every element a **named, recorded, swappable** choice with its
source, not a hardcoded default — because a head-to-head of tPA, ubiquitin and LAMP-1 found all three
beat the untagged vaccine while **none steered the arm it was chosen for** (PMID 19356616). Ship them
as options with that caveat attached, and let `deslip` run over the assembled CDS rather than the
cassette alone.

---

## Deliberately not on this roadmap

- **A processing predictor of our own.** NetChop/TAP/stability additions have repeatedly failed to
  improve on presentation prediction; V2.1 uses a published TAP motif, not a new model.
- **Nesting geometry** — see V3.3.
- **Duplication as a default** — see V3.4.
- **A better path search for the junction objective** — the free set is large and undifferentiated;
  V4.1 spends the freedom instead of narrowing it.

## The experiments that would actually settle this

Recorded because each is small, none has been run, and each collapses a `[silico]` or `[open]` row in
[`vector_evidence.md`](vector_evidence.md) into an `[exp]` one:

1. **`AAY` versus `AAA` in one processing assay.** The entire field uses `AAY` on a citation cascade;
   the one alanine result compared alanine-based against `GGGS`, not tyrosine against alanine. One
   SIINFEKL-style presentation assay decides it.
2. **A class-I/class-II mixed cassette with junction scanning both ways.** Nothing retrieved studies
   how the two junction types interact in one molecule, which is precisely what V1 builds.
3. **Copy number with a T-cell readout.** The non-monotonic copy-number result is B-cell display; the
   equivalent for a T-cell cassette has never been run.
4. **N versus 2N units at matched total dose.** This is the measurement that would fit `n0`, the
   module's one free parameter, and no trial has ever done it.
