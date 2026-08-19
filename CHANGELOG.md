# Changelog

All notable changes to `mhcmatch`. Format loosely follows [Keep a Changelog](https://keepachangelog.com);
versioning is [SemVer](https://semver.org).

> Note: 0.4.0–0.4.2 shipped without entries here. This file jumps 0.3.0 → 0.5.0; see `git log` for
> the 0.4.x range.

## [0.16.0] - 2026-08-19

### Added

- **`mhcmatch.complement` scores MHC class II.** It was class I only, because the class-I `aa` block
  bins its log-odds tables on the peptide's *length* and a class-II ligand is a 9-mer core floating
  inside an 11–25-mer. `score(peps, cls="mhc2")` takes its anchors from the P1/P4/P6/P9 core of the
  register (`store.anchor_indices`, or a frame you pin with `registers=`) and reads its own vendored
  tables, fitted on 603,781 human and 50,258 mouse class-II peptides. Class I is byte-identical and
  asserted so; the class is an argument and never inferred from the length.

  Which variable the block should be keyed on was measured rather than assumed, and the answer was
  not the predicted one — AUROC over the pooled role pair:

  | `aa` construction | human | mouse |
  |---|--:|--:|
  | register zones | +0.0029 | +0.0034 |
  | total length | +0.0070 | +0.0159 |
  | **both** | **+0.0102** | **+0.0185** |

  A class-II ligand's length is the length of its *flanks*, which is a covariate in its own right
  and not a register question. So both classes carry the same shape — pooled role pair, a length
  key, a position key — differing only in the position key: relative thirds of the TCR face at
  class I, register zones at class II.

- **`mhcmatch.vector.epitope_map` / `write_map` — the cassette, described.** One row per unit,
  linker and predicted epitope, 1-based over the cassette, as TSV and as JSON with a per-unit
  summary. Units and linkers tile exactly; an epitope spanning a junction is marked `unit = 0`; the
  class-II register core is resolved into cassette coordinates. **A peptide presented by two of the
  recipient's alleles gets two rows** — at a heterozygous locus those are two presentation events.
  `map_summary` reports per unit whether its class-I epitopes have overlapping class-II epitopes
  (`self_help`), the configuration Kissick et al. showed can replace an exogenous helper outright
  (PMID 24690990). CLI: `--map`, `--map-json`, `--map-threshold`, `--map-alleles-mhc2`.

- **`integrations/nextflow/mhcmatch/slurm.config`** — executor, per-process resources, scheduler-kill
  retries, and one shared `MHCMATCH_PMHC_DIR` / `MHCMATCH_CALIBRATION_CACHE` for the whole run.

- **`docs/safety.rst`** — the exclusion policy, the prior-evidence columns and what `n0` means, with
  the measurements behind each. Says two things that were nowhere written down: the prior-evidence
  columns are self-fulfilling on our own corpora and informative only on fresh data; and the screen
  is **class I / CD8 only** by design, CD4 self-reactivity being a different question.

### Changed

- **`complement`'s `motif` block is documented, and one claim in it was wrong.** A non-standard
  residue **breaks** a hydropathy run rather than being transparent to it — `AAAIIXIAA` gives
  `kd_run_max = 2`, exactly as `AAAIIDIAA` does — and it sits in `kd_run_frac`'s denominator without
  ever entering the numerator. The threshold is the median of the Kyte–Doolittle scale itself
  (−0.85, admitting `ACFGILMSTV`). The block's recorded gain is now in the docs: positive on all
  eight corpus arms, median +0.0060 AUROC.

- The Nextflow subworkflow runs `NEOAG`, `MIMICRY` and `VECTOR` on **class I only**, by design
  rather than omission; `PREDICT` and `RANK` still serve both.

- Container tag and conda pin moved from the 0.14.0 the module still named.

### Fixed

- **`fetch_pmhc` ignored the local mirror.** It called `hf_hub_download` directly, so
  `$MHCMATCH_PMHC_DIR` — the dataset root the SLURM profile and the cluster README export — was
  honoured by `fetch_file`/`fetch_proteome` but not by the one accessor every `Store.from_pmhc()`
  goes through. On a cluster following our own instructions each task therefore reached HuggingFace
  from a compute node instead of reading the staged mirror. It now resolves through `fetch_file`,
  exactly as `fetch_proteome` already did.
- **`Store.from_pmhc()` did not expand `~`.** A path like `~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz`
  raised `FileNotFoundError` on a file that was present — including the example in `skills/mhcmatch/SKILL.md`.
- **`slurm.config` declared its params after using them.** `process.queue` and the `env` block are
  plain assignments, so `queue`, `MHCMATCH_PMHC_DIR` and `MHCMATCH_CALIBRATION_CACHE` read back
  null: tasks went to the default partition and ignored the shared reference and calibration
  directories, silently.
- **Version pins were a release behind** in `main.nf`, `Dockerfile`, `environment.yml`, the module
  README and the `__init__.py` fallback (all `0.15.0`). Two tests now assert they track
  `pyproject.toml`, since this had drifted twice.
- **Docs corrected**: `complementarity.rst` still stated that no fitted class-II recognition model
  exists, which this release makes false — the warning is now scoped to `recognition.score_mhc2`,
  which remains MHC-I coefficients on a class-II core, and the fitted `complement.score(cls="mhc2")`
  has its own section. `api.rst` and `README.md` described the `aa` block as class-I only; the
  notebook count said six for seven.

## [0.15.0] - 2026-08-19

### Added

- **`mhcmatch.recognition` — the recognition head, as three models rather than one.** Each is
  fitted alone so their fit criteria are comparable and each score is readable on its own terms.
  The default is whichever wins BIC, currently `posbayes` for both species.

  | head | parameters | what it is |
  |---|--:|---|
  | `posbayes` | 3 | naive Bayes over amino-acid identity conditioned on **face**, scored as a summed log-likelihood ratio. Two 20-cell tables. Pure numpy |
  | `physchem_glm` | 23 | raw Kidera sums per face; `KF0` is the constant 1, so its face sums are the face sizes and length is never a separate feature |
  | `esm64_glm` | 65 | 64 components of a whole-peptide ESM2 pool. Most accurate on mouse, least explainable |

  **The default head needs no optional dependency.** A user who never installs `mhcmatch[esm]` gets
  a complete fitted model, not a degraded one.

- The split is by **face**, not by absolute position, because peptide length is not fixed and a
  model conditioned on position is not well defined across an 8-mer and an 11-mer.
- `recognition.log_odds_table()` prints the whole of `posbayes` — forty numbers.
- `score_mhc2` applies the MHC-I coefficients to the class-II binding core with P1/P4/P6/P9 as
  groove-facing. There is no fitted class-II model and no corpus to fit one on; it warns at runtime
  and the docs say so in a box. Scoring the core rather than the whole peptide keeps the face sizes
  inside the fitted range.
- `complement.kidera_design(peptides, anchors=…, roles=…)` — all ten Kidera factors by role.
- `tools/build_recognition.py`, carrying its own PEP 723 environment.
- Optional extra `mhcmatch[esm]`; the `esm64_glm` head raises if it is absent rather than dropping
  its features silently.
- `bootstrap --reference` fetches `immunogenicity/chowell_iedb_full.tsv.gz`.
- A `PROVENANCE.md` entry for `complement_mhc1_*.json`, which had none.

### Notes

- Coefficients come from `chowell_iedb_full_matched` — the rebuilt corpus with negatives resampled
  so the allele group carries no signal about the label. Measured cost against the unmatched arm,
  stated once: about 0.02 (human) and 0.06 (mouse) held-out AUROC.
- `mhcmatch.complement` is unchanged and still shipped. `recognition` is an addition; the recorded
  AUROCs for `complement` still belong to the arms it was fitted on.
- On the matched arm the two `posbayes` face tables correlate +0.94 (human) and +0.86 (mouse), with
  3/20 and 2/20 residues differing in sign. The face split is what makes the model length-agnostic,
  but on this corpus the two faces largely agree.

## [0.14.0] - 2026-08-18

The cassette gets a nucleotide half, the pipeline gets the rest of the library.

### Added

- **`vector.back_translate`** — the coding sequence for a cassette. Highest-usage human codon per
  residue from :data:`vector.CODON_USAGE_HUMAN` (Kazusa, *Homo sapiens* [gbpri], 93,487 CDSs /
  40,662,582 codons), backing off to the next synonymous codon whenever the first would extend a
  homopolymer past :data:`vector.MAX_HOMOPOLYMER`, then `deslip`. `mhcmatch vector --fasta-nt`.

  This is **not** a codon optimiser and does not claim to be. It fixes the two things that make a
  *polyepitope* fail where a natural ORF would not — the m1Ψ +1-frameshift motif, which a concatemer
  hits far more often because the designer chooses the seam residues, and synthesis-hostile
  homopolymers, which spacers like `AAA` manufacture directly. GC content, secondary structure,
  splice sites and CpG are untouched; a manufacturer's own optimiser should be preferred where there
  is one. The backoff is greedy, so `max_run` is a **target rather than a bound**: measured over
  5,000 random 20–60mers, longest run 6 and 84% at or below 4, against 13 for
  most-frequent-codon alone. Poly-proline pins the floor — all four proline codons begin `CC`, so
  consecutive prolines cannot be brought below a 5-run by any synonymous choice.
- **`vector.translate`** — so "synonymous" is checkable rather than asserted, by `deslip`,
  `back_translate` and by a caller supplying their own table.
- **`vector.units_from_context`** and **`mhcmatch vector --context windows.fasta`** — the join from
  `rank` to a unit table, which was the open item that made cassette assembly a manual step. `rank`
  emits minimal epitopes and a unit is the long window around the mutation; where that mutation sits
  lives in the FASTA header, so neither side alone can build one. Rows are grouped by **variant, not
  by register** — twenty registers of one mutation are one thing to put in a cassette, and `select`
  spends capacity per unit.
- **`rank.BASE_COLUMNS` / `EXTENDED_COLUMNS` / `ANNOTATE_COLUMNS` / `columns()`** — the `mhcmatch
  rank` schema as data, so a consumer can name the columns without running the command. Hoisted out
  of the CLI because a schema typed a second time is a schema that drifts.

### Changed

- **The nextflow integration covers the library, not just `predict`.** Four new processes —
  `MHCMATCH_RANK`, `MHCMATCH_NEOAG`, `MHCMATCH_MIMICRY`, `MHCMATCH_VECTOR` — plus
  `subworkflows/mhcmatch.nf` chaining all five, and a README written around each one's input and
  output contract. The image now bakes `bootstrap --reference`, because `rank`/`neoag`/`mimicry`
  read the known-epitope sets, mimicry references and expression tables and would otherwise reach
  for HuggingFace from a compute node.

### Fixed

- **The nextflow stubs emitted the wrong schema** — an 18-column `scored.csv` header against the
  real 57 and a 5-column `native.tsv` header against the real 27, so `-stub-run` produced files that
  did not match the ones a real run makes. No stub types a header any more: each asks the installed
  library for its own, and cannot drift again.
- Version pins in the nextflow module (`main.nf`, `Dockerfile`, `environment.yml`, README) were
  three minors behind at 0.10.0.
- The source-tree `__version__` fallback said `0.12.0`. It is what every `versions.yml` in the
  nextflow module reports, so a stale value mislabels a pipeline run.

## [0.13.0] - 2026-08-18

The step after ranking: assembling a cassette, and refusing one.

### Added

- **`mhcmatch.vector`** — polyepitope cassette assembly. `screen` **excludes** candidates on
  essential-tissue risk and runs *before* `select`, because capacity spent on a unit that has to be
  withdrawn is capacity not spent on a safe one. `select` grows each allotype while the next
  candidate beats that allotype's own expected yield per slot, so diversification falls out of the
  arithmetic rather than being imposed as a quota; `order` picks a spacer and an ordering minimising
  the strongest predicted binder spanning each junction, trying **no spacer first**;
  `slippery_sites`/`deslip` find and synonymously remove the m1Ψ +1-ribosomal-frameshift motif, which
  matters more for a concatemer than for a natural ORF. Scoring is injected (`binder`, `risk`), so
  the layout and policy logic are testable with no panel, no proteome and no download.
- **`vector.self_origin_risk`** — the shipped exclusion policy, two clauses: the unit's own target
  gene transcribed in an essential tissue (the MAGE-A12 shape), and a register **exactly** coinciding
  with an **unrelated** protein that is (the titin shape). Hits to the unit's own parent are excluded
  — a 27-mer is native context by design, and without that exclusion the screen rejects every unit of
  every cassette. Two defaults are measured rather than assumed:
  - `min_tpm=0.25`, not the conventional 5, because the two fatal precedents differ by two orders of
    magnitude and the lower one is what has to be caught — titin 64.4 TPM in heart left ventricle,
    MAGE-A12 **0.33** in brain caudate.
  - `max_subs=0`, because the decision is per unit while the search is per register and a 27-mer
    carries ~70 of them. At radius 1, 1 of 6 hazard-free random 27-mers is withdrawn at 9-mers and
    **4 of 6** across 8–11mers; radius 0 is clean at every length and still catches the titin unit
    (`bench/results/vector_screen_radius.md`).

  A mimicry-similarity screen was built and measured first and fires on almost everything — FPR
  0.693 against 0.020 at equal sensitivity — because anchor-masked similarity to a presented
  reference is presentation rather than recognition (`bench/results/vector_safety_screen.md`).
- **`mhcmatch vector`** — the whole cassette pipeline in one call:
  `--candidates units.tsv --n0 8 [--screen] [--fasta out.fa]`. Its input is a table of **long
  windows**, not `rank`'s minimal epitopes, and the reader says so when a column is missing rather
  than silently building the tolerising configuration. `--n0` is required with no default, matching
  the library: nothing in the public record fits per-allotype capacity. `--screen` is opt-in because
  it costs a whole-proteome index per register length — **without it no safety check runs at all.**
  The report names every withdrawal, every allotype's spend, every unselected candidate with the
  threshold it missed, and every junction.
- **`mhcmatch deslip`** — the m1Ψ +1-frameshift scan as its own command, since it takes nucleotides
  rather than peptides. `--fix` writes the repaired CDS, `TTT` → `TTC` upstream, protein unchanged.
  On a clean sequence it says so *and* says the check only applies to an m1Ψ construct.
- **`proteome.gene_symbols`** — `{name|accession: gene}` from the UniProt `GN=` field. Closes the
  join between deposits naming sources as accessions and `expression.safety_profile`, which is keyed
  on HGNC symbols. `mimicry.safety` takes it as `symbols=` and now returns a resolved `gene` beside
  the deposit's own `source`.

### Fixed

- **`expression.safety_profile` scanned all 5,586,792 rows per call** — 511 ms, and its callers ask
  per gene inside a loop (`mimicry.safety` once per mimic hit). Indexed by gene: 0.1 µs, identical
  values. The index is keyed on the resolved file rather than the `path` argument, so it cannot serve
  a stale table after `$MHCMATCH_EXPRESSION` changes.

## [0.12.0] - 2026-08-17

Mimicry stops being a distance and becomes a signed, per-component risk; the precursor estimators
move to the repertoire library that owns them.

### Added

- **`mhcmatch.mimicry`** — the fitted aggregate. Three references (`viral` priming, `self` tolerance
  *and* the autoimmunity read-out, `thymus` negative selection), each split into an **anchor** and a
  **TCR-facing** channel that partition the peptide, so no position is weighted twice. Six signed
  log-odds contributions and their sum, from a Bayesian logistic fit over 337,972 rows / 1,719
  positives across seven screens with screen indicators as nuisance columns — which is what makes
  the shipped coefficients within-screen. `mimicry_mhc1.json` v0.12.0.
  - **The signs follow the reference, as designed**: `viral` positive on both channels (+0.605
    anchor z = +16.8, +0.443 tcr z = +5.6), `self` negative on both (−0.304, −0.464), `thymus`
    positive on its anchor (+0.368) and unresolved on its TCR channel (+0.075, |z| = 1.1).
  - **A single whole-peptide distance was the wrong feature, and that was a search artifact, not
    biology.** Whole-peptide radius-2 thymic coverage of a candidate set is 1.63 % (viral 1.10 %);
    restricting to the TCR face at radius 1 reaches **53.4 %** against 0.25 % for the whole peptide
    at the same radius.
  - **Scores are log-odds; `probability()` is a separate step that demands a *named* corpus.** The
    seven screens run from 0.048 % to 46.8 % positive, so an unqualified probability mostly reports
    which intercept was used. AUROC **0.849 pooled / 0.596 median within screen** — both are in the
    artifact's fit record, and the second is the one to quote.
  - **Not collinear with the presentation stack**: max |r| 0.19 against affinity, 0.068 against
    agretopicity, 0.034 against expression, all VIF < 3.3. The TCR channel does correlate with
    `ipred` (r = 0.73–0.82), which is why its sign moves when `ipred` is already in the model — the
    module docstring keeps the two conditionings apart deliberately.
  - **`MimicryScore.nearest` carries which reference peptide was hit and what protein it came from**,
    so `mimicry.safety()` can resolve a self or thymic mimic through
    `expression.safety_profile`. Without the identity those channels are a bare number and the
    question a vaccine actually needs answered is unreachable. Thymic sources are UniProt
    accessions; `safety()` returns the accession with an empty profile rather than guessing a gene.
- **`mhcmatch mimicry`** — the aggregate on the command line, one column per (component, channel)
  plus the nearest hit and its source, so a rank can always be taken apart. `--corpus` for a named
  probability, `--no-self` to skip the expensive reference, `--coefficients` to print the shipped
  model and its fit record without scoring anything.
- **`mhcmatch neoag`** — annotate candidates against the tested-neoantigen database: nearest
  validated-immunogenic peptide and its substitution distance. With `--peptides` on a TSV every
  original column is carried through, so it drops into an existing candidate table without a join.
  **Prior evidence, never a fitted term** — every labelled screen we hold sits inside that database
  (retrieval recall at distance 0 is 1.000 on all seven), so a coefficient on it would be
  memorisation. Held out honestly it still earns its place: rebuilt without the test screen, fuzzy
  matching at two substitutions recovers 0.08–0.34 of a screen's positives where exact lookup
  recovers 0.00–0.26, which is why `--max-subs` defaults to 2 rather than 0.
- **`mhcmatch rank --extended` and `--annotate`** — the mimicry read-out on the candidate table
  that a donor's shortlist actually comes out of. `--extended` appends the six signed contributions,
  their sum and the autoimmunity total; `--annotate` appends what each candidate *resembles* — the
  nearest self / viral / thymic mimic per channel with its source protein, plus the nearest
  validated neoantigen and its distance. **Both are columns, never a re-score**: the base schema is
  a strict prefix and the ordering is byte-identical with and without them, asserted in the test
  suite, because whether mimicry belongs inside the gate is a benchmark question that is not
  settled. `--no-self` trades the expensive reference for speed on both.
- **`notebooks/07_mimicry_risk.py`** — the fitted form, against notebook 4's raw scan.

### Changed

- **`mhcmatch.precursor` is now a re-export of `vdjmatch.precursor`** (801 lines → 75). The
  estimators, their maths and the `vdjmatch precursor` CLI live in the repertoire library, which is
  where the problem belongs; this name keeps working so existing imports and notebooks do not break.
  Two behaviour changes came with the move, both because nothing enumerates any more: `shell_profile`
  has no `max_members` and no memory ceiling (the `r=2` profile for 300 junctions cost ~9.9 M
  materialised strings and is now a few DP passes), and its shells report `n=None` when the union is
  too large to census — the masses are always exact. Adds `union_mass`, `closed_ball_mass`,
  `unseen_junctions` and `precursor_frequency`. Install unchanged: `pip install 'mhcmatch[precursor]'`.

### Fixed

- **`mimicry.score()` silently produced a different, smaller model** when a component was missing
  from `refs` — the absent feature standardizes to zero rather than erroring, and the usual way to
  get there is `load_references(with_self=False)`, which drops the component carrying the largest
  coefficients. It now raises, with `allow_missing=True` to accept that deliberately.
- **A bare `|r|` in the `mimicry` module docstring** was an undefined RST substitution and failed the
  `-W` docs build.
- **The precursor tests guarded on `vdjtools` while the import is now `vdjmatch`** — backwards for
  the one case that matters, since anyone upgrading from 0.11.0 with the `[precursor]` extra has
  vdjtools already and vdjmatch not yet, so the guard passed and the import then failed. Guarded on
  the module actually imported.
- **The vendored anchor models are regenerated under 0.12.0.** They are version-stamped and a stale
  one is ignored and refitted at runtime (~200 s); the staleness test exists to fail a version bump
  that forgets them, and on this bump it did.

## [0.11.0] - 2026-08-17

A CLI that can be pointed at a file, a length-aware recognition model, and mimic categories that say
what they mean.

### Added

- **`--peptides FILE` on every peptide-keyed command** — `decompose`, `restriction`, `affinity`,
  `binder`, `source`, `explain`, `complement`, and the new `mimics`. One peptide per line or a TSV
  with a `peptide` column; `-` reads stdin; output is TSV on stdout or `--out`. This is the fix for
  the real cost of the CLI, which was never the scoring: the presentation and affinity calibrators
  are ~5 s, the binder calibrator ~45 s and a human-proteome length index ~70 s — all cached for the
  life of the process, and all re-paid by a shell loop. `bench/cli/run_cli_bench.zsh` measures both
  forms of every command.
- **`--threads` on `source` and `mimics`**, and deliberately nowhere else. Those two run the
  neighbour search in C++ with the GIL released; every other command's per-peptide work is a small
  numpy product, so the flag is absent rather than accepted and ignored.
- **`mhcmatch mimics`** — the module had no CLI at all. Reports near-identical reference peptides
  per category with that category's *kind*, batched and threaded.
- **`mhcmatch.proteome.Proteome.find_sources`** — the batch form of `find_source`: one index build
  per length and one threaded `search_batch`, instead of one Python-level query per peptide. Also
  `windows(L)`, the public form of the window set the mimic loaders need.
- **`mhcmatch.mimics.KINDS` and `PROTEOME_REFS`** — the mimic categories are now `thymus`, `self`,
  `viral`, `bacterial` and `neoag`, each with what a hit in it *argues*, and they are never summed.
  `self` (the host proteome) is kept separate from `thymus` (the thymic immunopeptidome) because
  being encoded does not imply being presented and the two license different conclusions.
  `bacterial` is five gut-commensal and pathogen proteomes; `load_reference_sets(..., proteomes=…)`
  builds them — class I only, since class II spans 15 lengths and would materialise tens of millions
  of windows, which now raises rather than swapping.
- **`affinity --peptides` reads a `wt_peptide` column**, so agretopicity comes out of the same pass
  instead of a second run joined back on a peptide string that is not a key.
- **`expression.TUMOR_TISSUE` and `matched_tissues()`** — each of the 19 tumour types now maps to its
  matched normal tissue, and `expression --list-contexts` prints the pairing instead of two unrelated
  lists. The safety read previously required the caller to already know that melanoma pairs with
  skin. Both vocabularies are named in the docs because **neither is clinical**: `--tumor` takes
  **TCGA study abbreviations** (NCI GDC), `--tissue` takes **GTEx `SMTSD`** names, so a pipeline
  needing ICD-O-3, SNOMED CT or OncoTree must bring its own crosswalk. `CRC` is flagged as not a TCGA
  code (TCGA has `COAD` and `READ` separately; the source table merged them) and `HNSC` as
  approximate (GTEx has no head-and-neck mucosa). A test asserts all 19 resolve, so a tissue name
  that rots fails loudly rather than silently emptying the safety read.
- **`pseudoseq.class2_report(key, mode)` and `--mhc2-report {pair,beta,isotype}`** on the five
  commands that *choose* an allele (`restriction`, `binder`, `scan`, `predict`, `rank`). A class-II
  key does not lead with the same chain at every isotype: DRA is monomorphic so DR is keyed by its
  **beta** (`DRB1_0101`), while DP and DQ are keyed by the alpha–beta pair and lead with the
  **alpha** (`HLA-DQA10501-DQB10301`). Anything that compares two callers by the leading gene is
  therefore matching DR's beta against DP/DQ's alpha, and splitting DR against itself whenever the
  DRB gene differs. Measured on a 10,402-row class-II concordance: leading-gene agreement 0.401,
  true DR/DP/DQ agreement **0.527**, the gap being 1,318 DR-vs-DR pairs. `pair` is the default and
  its output is unchanged; commands *handed* an allele still echo what the caller typed, and class-I
  and mouse keys come back untouched rather than reduced to a stub.

### Changed

- **The `aa` block of `complement` is length-aware**, and both species are refitted. It keeps the
  pooled `aa_anchor`/`aa_tcr` pair — whose sum is still exactly `posbayes.llr` — and adds an
  anchor/TCR table **per length bin (8, 9, 10, 11+)** plus the TCR face in **relative thirds**.
  Bins rather than one table per observed length, so a 12- or 13-mer is scorable at all. Against the
  pooled construction under peptide-grouped CV, paired bootstrap over peptide groups, CI excluding
  zero on all four corpus arms: chowell/human +0.0069, chowell/mouse +0.0115, kesmir/human +0.0206,
  kesmir/mouse +0.0208 AUROC. A length × role interaction and a bulge/flank split both buy nothing,
  which localises the effect: length carries *which residue is preferred where*, not a global
  reweighting. 19 features → 30. `bench/results/length_roles.md`.
- **`rank.GATE` refitted** on the new recognition axis — `recog_mu`/`recog_sd` describe *that* axis,
  so they move with it. Holdouts unchanged within noise: TESLA 0.592, Neopep 0.804, Gfeller 0.784.
- **Nextflow module pins bumped 0.8.0 → 0.10.0**, and the Dockerfile now states what its `bootstrap`
  does *not* cover, so an offline `rank` process fails at build time rather than on a compute node.
- **README rewritten** around a task → command table, the batch/threads contract, and the two axes.

### Performance

- **The binder calibrator was 45 s and is now ~4 s**, from two independent causes the CLI benchmark's
  profile exposed. Nothing it reports changes: the isotonic step levels are identical on 300
  randomized trials, and the Potts scores are **bit-identical** on 60,000 checks over three alleles.
  - `calibrate._isotonic` was **O(n²)**. PAVA is linear, but the blocks lived in a list and each pool
    did `del ys[i + 1]`, shifting the tail of three lists. A common allele's known ligands against a
    10,000-peptide background is ~118,000 points: 2.66 s → 0.041 s. Blocks go on a stack now.
  - `PottsAffinity.predict_y` summed ~315 weights per peptide, and ~34 of every 35 depend only on
    the allele. The pocket side is fixed once the allele is, so the energy factors into a constant
    plus a table `E[p][r]`, leaving nine float adds per peptide — **21× scalar, 30× via the new
    `predict_y_batch`**. Each cell is a `math.fsum`, and that is load-bearing rather than
    fastidious: the weights are float32 and the loop it replaces added them into a Python float,
    i.e. exactly. Summing the pocket contributions in float32 costs ~1e-7 and moves 735 of 20,000
    IC50 values at their reported precision; numpy's float64 pairwise sum leaves ~2e-9 and moves 122.
  - End to end: `mhcmatch binder` on one peptide 52.7 s → 10.7 s, `explain` 37.9 → 7.7,
    `predict` 48.2 → 9.5, `rank fasta` 53.3 → 9.6; the library's own test suite 71 s → 26 s.
    Commands that do not touch the calibrator are unchanged at 1.0–1.1×.

### Fixed

- **`mimics.DEFAULT_REFS["neoag"]` pointed at `immunogenicity/neoag_tested.tsv.gz`**, which 404s —
  the deposit moved it to `neoantigens/`. The documented default reference set was unusable without
  a local mirror. A regression test now asserts every default path resolves.
- **`store.fetch_proteome` ignored `$MHCMATCH_PMHC_DIR`**, so a local mirror was bypassed for
  proteomes only. Routed through `fetch_file` like everything else.

## [0.10.0] - 2026-08-17

The recognition axis, grown up: a six-block complementarity score with a per-species table, a
neoantigen ranker that no longer applies its coefficients to the wrong scale, built-in known-epitope
lookup, and a batched mimic search that is three orders of magnitude faster than the one it replaces.

### Added

- **`mhcmatch.complement`** — complementarity, i.e. how well a presented peptide complements a TCR
  repertoire. Six feature blocks: `ipred`'s PC1/PC2 and length; the same components split
  **MHC-facing vs TCR-facing**; **MJ1996** on the anchors and **TCRen marginalised over 28,250,990
  real TRB CDR3 loops** on the TCR face; contiguous-hydrophobic-run motifs; per-role **residue
  log-odds**; adjacent TCR-facing dipeptides. Emits a prior-free log-odds, with `posterior()` for a
  probability at the caller's own base rate.
  - The `aa` block's two columns sum to `posbayes.llr` **exactly** (asserted in the test suite), so
    that model is a strict special case and the block ablation measures what the other five add.
  - Beats it on all four deposited corpus arms × both hosts (chowell/human 0.7125 vs 0.7111,
    chowell/mouse 0.7633 vs 0.7582, kesmir/human 0.6480 vs 0.6369).
  - **Per species, never pooled**: `score(peps, species="mouse")` uses the 47,140-row mouse arm.
    Different MHC, different thymic repertoires — one fit across them is fitting a mixture.
  - The head is linear because a diagonal-covariance Gaussian **cannot represent a summed
    log-odds**; the EM and supervised Gaussian parameters ship alongside so the comparison stays
    re-checkable.
  - Vectorised: **511,301 peptides in 0.93 s**. The dipeptide block is a sparse `(code, row)` list,
    not a dense `(n, 400)` matrix.
- **`mhcmatch.known`** — five built-in reference sets for exact-match lookup, assembled from the
  public deposits: `neoantigen` (23,299 confirmed immunogenic tumour neoantigens from NCI/Gartner,
  the epitope-resolution screens and the aggregated cohorts), `neoantigen_neg` (468,220 screened and
  found non-immunogenic — the one label that says this exact peptide was tried and did not work),
  `immunogenic` (15,889), `self` (53,878 thymic), `viral` (44,993). `rank` uses them by default.
- **`mhcmatch.mimics.neighbours`** — batched same-length mimic search, and `scan(evalue=False)` to
  route through it. **237,000 queries/s against 55** for the per-query `find_mimics` path it
  replaces, on measured identical counts and distances.
- **`mhcmatch complement`** CLI — scores peptides or a whole TSV, `--features` to take a score
  apart, `--prior` for a probability, `--species`.
- **`store.fetch_file`** — any file of the public dataset by repo-relative path, so a worked example
  can run on a whole published deposit; `bootstrap --reference` pre-stages all six in one call.
- `docs/complementarity.rst`, `notebooks/06_complementarity.py`, `tests/test_complement.py`.

### Fixed

- **`rank.GATE` applied z-score coefficients to raw axes.** The fitting script standardizes both
  axes and never wrote the standardizer out, so `GATE` carried `mu = 0, sd = 1` placeholders. A
  product of two sigmoids is **not** rank-preserving under a monotone rescaling of one axis, so this
  moved the ranking and not merely the calibration. Refitted with the standardizer recorded: every
  cohort improves — TESLA 0.597 vs 0.473, Neopep 0.802 vs 0.662, Gfeller 0.782 vs 0.702 AUROC.
- **`mimics` `n_near` counted deposit rows, not peptides.** The compendia repeat a peptide once per
  allele/source it was reported under (the viral set is 57,331 rows over 26,640 distinct), so the
  count was a function of deposit frequency rather than of the sequence neighbourhood. The batch
  path deduplicates; `top_mimic` and `top_subs` are unaffected.

### Changed

- **`rank`'s recognition axis is now `complement.score`**, not `posbayes.llr` — measured on
  peptide-grouped CV over every corpus arm and host. `mhcmatch explain` prints both, plus `ipred`.

## [0.9.0] - 2026-08-16

Three new public modules on the recognition side of the problem — physicochemical featurization, a
frozen immunogenicity model, and TCR precursor frequency — plus a calibrated probability on the
binder score. Nothing on the presentation path changes behaviour; the vendored anchor models are
regenerated only so their version stamp matches this release (panel unchanged).

### Added

- **`mhcmatch.immuno`** — physicochemical featurization of an epitope: 141 features per peptide
  (20 amino-acid scales × 7 statistics, plus `length`), with the two contested choices exposed as
  arguments rather than baked in. `ANCHOR_SCHEMES` keeps all three class-I anchor definitions in the
  toolchain selectable, and the `sum`/`mean`/`min`/`max` descriptors are joined by
  `run_max`/`run_n`/`run_frac`, which express *contiguity* — a property no composition statistic can
  represent. Needs no reference panel and no download. Self-check: `python -m mhcmatch.immuno`.
- **`mhcmatch.ipred`** — the fitted physicochemical immunogenicity model over that basis: two
  principal components of the property matrix plus length, thirteen parameters, returning a
  calibrated `log P(immunogenic)` (`ipred.p_immunogenic`). Parameters are vendored in
  `mhcmatch/data/ipred_mhc1.json` and never refitted at import time.
- **`mhcmatch.precursor`** — TCR precursor frequency `F(e)` for an epitope: six independent
  estimators (`observed_mass`, `coverage_corrected_mass`, `ball_mass`, `shell_profile`,
  `event_ratio`, `motif_mass`) plus `cross_check`, which is the point — they bound the answer from
  different directions instead of agreeing by construction. Optional extra,
  `pip install 'mhcmatch[precursor]'` (needs `vdjtools`); `check_junctions` guards the
  CDR3-vs-junction trap before any Pgen is computed.
- **`BinderScore.p_binder`** — a calibrated `P(binder)` alongside `binder_rank`. The %rank is what
  you sort by; `p_binder` is what you threshold or hand to a downstream model, because it means the
  same thing outside the candidate list it was computed in. Isotonic-fit from the allele's own
  ligands against the random-peptide background when the calibrator is built with `positives=`.
- **Structure-derived contact profile** (`mhcmatch.data.contact_profile`, reached as
  `immuno.contact_profile` / `scheme="contact"`) — continuous per-position TCR-facing weights from
  8,062 TCR↔peptide residue contacts over 370 crystal structures, with both derived steps (zeroing
  below half the uniform-footprint expectation, rescaling survivors to mean 1) fixed by the profile
  rather than tuned. On class-I 9-mers it recovers P1/P2/P3/PΩ as anchors unsupervised — which is
  neither shipped anchor scheme.
- **Four marimo notebooks** (`notebooks/`, `pip install 'mhcmatch[notebooks]'`) — presentation and
  the binder score, immunogenicity features, precursor frequency, mimicry and self. Clone-only (the
  wheel ships `src/mhcmatch` alone), but each bootstraps its data from HuggingFace or from the
  vendored tables, so none needs a local file.

### Documentation

- **Immunogenicity features** ([`docs/immunogenicity.rst`](docs/immunogenicity.rst)) — install to
  feature matrix, with the four position schemes side by side. The featurizer previously appeared
  nowhere on the docs site.
- **The amino-acid property basis** ([`docs/property_basis.rst`](docs/property_basis.rst)) — two
  properties of the vendored tables in `mhcmatch.data.aa_tables`, each pinned by a regression test:
  the dominant eigenvector of the 20 × 142 residue-by-scale matrix is a hydropathy axis (32.79 % of
  the variance; median |ρ| 0.894 against 39 named hydrophobicity scales), and the Kidera factors are
  already orthogonal (largest off-diagonal correlation 0.0026, participation ratio 10.00 of 10), so
  PCA over the alphabet returns an arbitrary rotation and reduces nothing. Scoped deliberately: that
  degeneracy holds under the uniform measure over residue types and breaks under any other.
- `docs/api.rst` gained the five modules it was missing — `immuno`, `ipred`, `precursor`,
  `data.aa_tables`, `data.contact_profile`.

## [0.8.0] - 2026-07-18

Gamaleya/ISPRAS beta-test feedback (170726), plus the generalized binder score.

### Added

- **Generalized binder score** (`store.binder_score` / `mhcmatch binder` / `predict.binder_score`) — a
  **calibrated combined %rank** fusing the presentation %rank (`AnchorModel`) and the affinity %rank
  (`PottsAffinity`): Fisher's combined statistic `-(ln p_pres + ln p_aff)`, itself calibrated per allele
  against a random-peptide background so `binder_rank` is a true %rank (correctly banded, cross-allele
  comparable). A soft-AND — scores well only when a peptide is *both* presented and binds. The two heads
  disagree along the binding-strength axis (presentation rescues weak-but-presented ligands, affinity
  rescues strong-but-atypical binders; Spearman(Δ, log nM)≈+0.5–0.65 on TESLA/NCI), so the blend is more
  robust than either alone — combined immunogenicity AUROC beats both single heads (TESLA 0.786, NCI 0.965).
- **Binder score flows through the pipeline.** `predict_windows` now annotates every predicted binder
  with `affinity_rank`, `binder_rank`, and `binder_band`, and `write_native` emits them — so the
  Nextflow module's `.mhcmatch.native.tsv` carries the generalized binder score with no extra call
  (fixed ~10 s one-time calibrator fill, cached per store). The `.scored.csv` keeps its fixed 57-column
  pipeline schema untouched.

### Fixed

- **`setup.sh` was fish-only.** Rewritten in POSIX shell so it runs under **bash, zsh, or sh**
  (calls `.venv/bin/pip` directly; no `source …activate`); `README.md` and `docs/getting-started.rst`
  invoke it as `bash setup.sh` again.
- **Quickstart referenced a non-shipped file.** `Store.from_pmhc("pmhc_full.tsv.gz", …)` →
  `Store.from_pmhc(tier="shortlist", …)` (auto-fetched from HF). `from_pmhc` now raises an actionable
  `FileNotFoundError` (pointing at `tier=` / `$MHCMATCH_PMHC`) instead of a bare `open()` error.
- **`StructureScorer` hard-coded a personal template path.** The default template dir was a fixed
  `~/vcs/code/tcren-ms/data/Canonical2026`, so a missing `1oga.pdb.gz` broke it. It now resolves via
  `tcren`'s own `data_dir()` (`$TCREN_DATA_DIR` or an editable checkout), keeps the
  `$MHCMATCH_STRUCTURES` override, and raises a clear error when a template PDB is absent.
- **MHC-II `predict` on a large input "never finished."** The register + K=3 motif EM (~200 s on the
  full corpus, paid twice per run) is now shipped **pre-fit** in `mhcmatch.data` and loaded read-only
  by `Store.anchor_model`, guarded by version + panel hash + build params. Loaded models are
  bit-identical to a fresh build (no benchmark number changes); a 1034-window MHC-II sample now runs
  in ~27 s instead of never. Read-only vendoring avoids any cache race under concurrent (nextflow/
  SLURM) execution. Both classes are shipped so the version/panel-hash guarantee is uniform. The
  release workflow (`publish.yml`) **regenerates the models before building the wheel**, so a published
  release can never ship stale models; `ci.yml`'s staleness test is the earlier (data-free) guard.
  Regenerate manually with `python tools/build_anchor_models.py`.

## [0.7.2] — 2026-07-17

**Three global constants were wrong on a heterogeneous panel; two now have per-allele/per-position
estimators.** Every knob below **ships inert at its default and is measured byte-identical**, so no
committed number re-baselines and nothing is a behaviour change until it is opted into. The headline is
diagnostic rather than a default flip: the class-II frequent gap is a **register-EM convergence failure
on HLA-DP**, not a motif deficit or an estimator-variance problem.

Results: `2026-mhcmatch-benchmark/KEY_FINDINGS.md`, `bench/results/register_em_convergence_dp.md`,
`bench/results/blosum_pseudocount.md`. Design: [`docs/hierarchical_rules.md`](docs/hierarchical_rules.md).

### Added

- **`AnchorModel(register_em="converge")`** — run the best-frame register EM to convergence **per
  allele** (freeze each one when its own frame assignments stop moving) instead of a shared pass count.
  No count serves the panel: HLA-DP is still improving at 32 passes while the rare stratum reaches its
  fixed point by 8 and never moves again, so the shipped `2` is an *early stop that flatters rare*, not
  a correct value. Measured on MHC-II human screening (K=3): frequent AUPRC **0.625 → 0.667**, gap to
  NetMHCIIpan-4.3i **−0.149 → −0.108 (28% closed)**. It **dominates every constant tried** — equal to
  `em=32` on frequent, better on medium (0.510) and rare (0.635), and **1.36× cheaper** (73 s vs 100 s),
  because frozen alleles skip the frame search.
  - The gain is DP-specific (**+0.043** mean vs DR **−0.005**) and the causal test passes:
    HLA-DPA1\*01:03/DPB1\*04:01, the DP allele already converged (H/Hmax 0.635), moves **+0.000 exactly**.
    No threshold, no allele family named, no benchmark label — DP earns its passes by still moving, and
    DRB1\*04:04 (0.2% eluted-ligand, boundaries genuinely arbitrary) keeps its flat prior rather than
    being forced to sharpen.
  - **Not the default, deliberately:** it is a *screening* win and a *restriction* cost — on
    `--decoy-mode hard` frequent barely moves (+0.001) while rare PPV@P flips from a win over
    NetMHCIIpan (0.402 vs 0.372) to a loss (0.350). A knob that must flip per task is usually still
    wrong; see `hierarchical_rules.md` for the frame-tally fix that should remove the trade.
- **`AnchorModel(prior_strength="auto")`** — empirical-Bayes shrinkage concentration **per anchor
  position**, by method of moments on the Dirichlet-multinomial
  (`τ_j = Σ_r m_j(r)(1−m_j(r)) / Var_between(j) − 1`), estimated on alleles with n ≥ 200 (where sampling
  noise is negligible) and applied to all. One global `τ=10` is wrong **in opposite directions at once**:
  between-allele PWM variance spans **71×** across MHC-I core positions, so at P4 (alleles barely differ)
  τ=10 leaves 33% of a rare allele's sampling noise in, while at P2 (alleles differ enormously) it
  discards 67% of its only real signal.
  - **Recovers the known anchors unsupervised**, which is the check that it measures what it claims:
    MHC-I P2 τ=**1.0** (B pocket) and PΩ τ=**1.7** (F pocket) against P4 τ=**71.5**; MHC-II's four lowest
    are P1/P4/P6/P9 — the hardcoded `MHC2_ANCHORS`. The global τ=10 is correct for **exactly one position
    in nine** (MHC-I P3). MHC-II's spread is 6× where MHC-I's is 71×: the open groove as a number.
  - Measured: MHC-II screening **rare AUPRC 0.648 → 0.689 (+0.041)**, extending the margin over
    NetMHCIIpan from +0.038 to **+0.079** (PPV 0.534 → 0.594) — the largest rare gain measured. It acts
    where τ carries mass (67–77% at rare, 0.9% at frequent). MHC-I restriction frequent holds and nudges
    up (AUPRC 0.850 → **0.854**); rare 0.749 → 0.726 flips to a loss.
  - **`converge` and `"auto"` do not compose**: together they keep the frequent gain (0.668, best PPV
    0.629) but τ's rare gain vanishes (0.689 → 0.630). That is a *positive* result about the mechanism —
    τ fixes **residue** borrowing while rare's damage under convergence is in the **frames**, which are
    tallied at full weight though the model that chose them was 67–77% borrowed. It locates the next fix.
  - Lengths and core offsets keep a scalar (`_tau_scalar`): they are not residue distributions, so a
    per-residue-position τ is meaningless for them. τ is fit on the **final** prefs, after the register
    EM (which bootstraps on the scalar), so the EM, the background null and the mixture assignments are
    unchanged.
- **`AnchorModel(pseudocount=β, pseudo_matrix=None)`** and **`pseudoseq.blosum62_conditional()`** — a
  mass-preserving BLOSUM62 substitution pseudocount on the anchor counters, `ĉ(r) = (1−w)·c(r) +
  w·Σ_r' c(r')·P(r|r')` with `w = β/(n+β)`. The Nielsen et al. 2004 recipe (PMID 14962912) that
  NetMHCpan's own lineage has used since 2004 and mhcmatch never had. **Ships off (β=0) because it is a
  measured negative** — see below. `P(a|b) = p_a·2^(s_ab/2)` needs no q_ij table and no new dependency
  (seqtree's BLOSUM62 was already imported for the allele kernel).

### Measured and rejected (recorded, not shipped)

- **BLOSUM pseudocounts make class-II screening monotonically worse**: frequent AUPRC 0.625 → 0.622 →
  0.618 → 0.612 → 0.602 over β = 0/25/50/100/200; the gap *widens* −0.149 → −0.173. The premise was sound
  and stands — only 28.0% of *frequent* MHC-II (allele, anchor) cells observe all 20 residues, and the
  count-0/count-1 boundary is a **3.8-nat cliff on a ~1σ Poisson difference** (HLA-A\*30:01 P2, n=734).
  **Mechanism, pre-registered before the run:** grading the never-seen penalty improves *bulk* ordering
  (rare/medium AUROC +0.006/+0.009 at β=25) but lifts the chemically plausible **near-miss** decoys that
  sit at the **top** of the ranking — which is what AUPRC and PPV measure. Every screening decoy is a
  proteome window, so its residues are plausible by construction. **The model's overconfidence about
  never-seen residues was doing useful work.** This ruled out estimator variance and redirected the
  search to the register.
- **MJ contact potentials not adopted**: measured **79% rank-1** (essentially a hydrophobicity axis), so
  they cannot express "an R pocket takes K but not S", and they need a temperature unsettable from first
  principles — where BLOSUM's conditional is parameter-free (reproduces the matrix to KL ≤ 0.011
  bits/column, argmax agreeing in all 20 columns; recovered `q_ab` symmetric to 5.1e-04). `pseudo_matrix`
  exists so the bench can pass an MJ conditional without mhcmatch vendoring MJ data or taking a `tcren`
  dep.
- **`eps=1e-3` is not the lever**: it *does* extinguish the τ prior at frequent alleles (prior mass
  1.25e-05, ~80× below eps) and clips decoys asymmetrically (13.7% of MHC-I frequent decoy lookups vs
  0.3% of positives) — but the metric is **flat from eps=0 to 1e-3**. Clipping shifts decoys roughly
  uniformly, and uniform shifts do not move a ranking. Left exactly where it is.

### Docs

- [`docs/hierarchical_rules.md`](docs/hierarchical_rules.md) — the design: global prior → family
  (kernel communities, Q=0.94/0.90) → allele, with the shrinkage strength derived from the variance ratio
  rather than tuned. Names the remaining violator: `footprint`'s `rare_max=30`, a capacity threshold
  sitting **exactly** on the evaluation stratum's boundary.
- `ROADMAP.md` §6b — the presentation-null item is **mostly shipped**, not open (`background="proteome"`
  is the `log(θ_A/p_proteome)` it prescribes, it is the CLI default, and the screening benchmark has been
  running it all along). Records the three refuted mechanisms so no future session re-chases them.

## [0.7.1] — 2026-07-17

**Potts affinity weights refit under the de-duplicated 8-mer encoding.** A correctness release: it
activates the `enc=1` fix that has been dead code since v0.6.1, and makes the vendored weights
reproducible from a documented command. **Every MHC-I and MHC-II affinity number changes.** It is
**not** a performance release — the refit is neutral within noise, measured, and that is on the record.

### Changed

- **`data/affinity_potts_mhc{1,2}.npz` refit** (`meta[4]=1`). MHC-I 22,971 → 29,651 nonzero weights,
  `b` +0.1185 → +0.0003; MHC-II 30,929 → 31,551, `b` +0.2819 → +0.1875. Two things move together and
  neither is a method change:
  - **The 8-mer collision is now actually fixed.** v0.6.1 fixed the *code* on both sides and bound the
    encoding to the weights via `meta[4]`, so the fix could only activate atomically with a refit —
    which never came. Every shipped 8-mer score until now used the legacy `core[:5] + core[-4:]` slice,
    where index 4 fills two slots and contributes two perfectly-correlated field terms. **8-mer scores
    change materially; L≥9 scores change only via the refit below.**
  - **The training set grew 73,880 → 84,709 points / 108 → 132 alleles.** The weights were fit
    2026-07-15 against `mhci_pseudo.fa` naming **4,143** alleles; `3bda000` ("68% of alleles were
    unscorable") and `0cd2d42` ("+7,085 alleles") landed **the next day** and took it to **20,082
    names / 5,407 grooves**, and the weights were never refit. All 4,143 old keys carry a
    byte-identical 34-mer today (0 changed, 15,939 added) — the fix *added* alleles, so the old weights
    were under-trained, never wrong.

  This also **resolves the "shipped weights are unreproducible" note** in the benchmark repo's
  `results/potts_mhc1_encoding_defects.md` (shipped 22,971 nonzero vs a fresh refit's 29,666 *with the
  legacy encoding restored*). The cause was the pseudosequence table, not `measured.tsv` drift; the old
  weights reproduce bit-exactly under `mhci_pseudo.fa@9e2444f`. Nothing needed pinning.

### Added

- **Regression tests for the vendored weights** (`tests/test_affinity.py`) — `meta[4] == 1` per class, an
  8-mer slot-mapping assertion, and pinned IC50 values for three (peptide, allele) pairs. There were
  **none**: a weight swap or a silent refit changed every shipped affinity score and still passed CI.

### Measured, and deliberately NOT shipped

- **BLOSUM/MJ "smarter than one-hot" encoding — tested, null, dropped.** `train_potts.set_soft(tau,k)`
  had implemented BLOSUM admixture on the groove axis all along, pinned to one-hot, never swept. Swept
  jointly with `alpha`, paired, 5 seeds: every arm lands inside **±0.010** rho against a 0.166
  common−rare gap. The reason is structural, not a shrug — soft encoding is *generalized ridge* under
  metric `(SSᵀ)⁻¹` (verified to 2.2e-16), and `S` is full-rank at every `(tau,k)`, so it adds **zero**
  new directions. Predicted to act like `alpha ×2.5`; measured, soft(τ=2,k=5)@α=40 reproduces
  one-hot@α=80 to within noise. `alpha=40` is already optimal, so there is nothing to win. Softening
  the *peptide* axis (which the design pins hard, and which NetMHCpan-4.0 does not) is the only arm with
  consistently positive signs and it is worth **+0.004**. Full result and mechanism:
  `bench/results/potts_encoding_ablation.md`.
- **Defect 1 (length-blindness) is still live and still unfixed.** `SLYNTGATL` and `SLYNTAAAGATL` score
  bit-identically. Per-length intercepts were measured here and are null on per-allele Spearman: the
  large effects (8-mers bind **5.5×** weaker than 9-mers within an allele) sit at 5.6% of the corpus.
  The recorded **+0.059 AUROC** for a length prior belongs to the *NCI immunogenicity ranking* task, not
  affinity regression. Tracked in ROADMAP §6c.

### Fixed

- `bench/affinity/fit_potts.py` wrote to `MultiplexedPath('…')` as a literal directory name when `--out`
  was omitted (`mhcmatch.data` is a namespace package, so `str(resources.files(...))` is a repr, not a
  path) — the default target never worked. *(benchmark repo)*

## [0.7.0] — 2026-07-17

**Per-allele motif mixtures for MHC-II, on by default.** A class-II allele now scores a mixture of
`K` PWM components (`AnchorModel(n_motifs=3)`, the new default) instead of one, closing ~40% of the
frequent-stratum AUPRC gap to NetMHCIIpan-4.3i. No API break — `n_motifs=1` restores the single-PWM
model and never enters the mixture path. MHC-I is unaffected (the mixture is class-II only).

This is the other half of GibbsCluster-style deconvolution: v0.6 marginalised over the binding
*register*; this fits the *motif*. It answers the "can extra matrices help?" question — and the
answer is a mixture, because the score is a sum of per-position log-odds and that family is closed
under addition, so any additive "extra matrix" collapses to one PWM. Only `log Σ_k π_k exp(s_k)` adds
capacity.

### Added

- **`AnchorModel(n_motifs=K)` / `Store.anchor_model(n_motifs=K)`** — K motif components per allele,
  fit by EM on the whole corpus (no external labels, no NetMHCpan), scored as
  `log Σ_k π_k Σ_r P(r|L,a)·exp(s_{k,r})`. Default **3** for MHC-II. Capacity self-adapts with **no
  ligand-count threshold**: a component with no counts for an allele returns that allele's pooled
  (shrunk) motif *identically*, so a thin allele degrades to the single PWM. Symmetry is broken by a
  deterministic `crc32(peptide) % K` init (reproducible; no seed to plumb).

### Changed

- **MHC-II scoring uses the K=3 mixture by default.** Measured, human MHC-II holdout (seed 0), frequent
  stratum AUPRC vs NetMHCIIpan-4.3i: allele-specificity **0.558 → 0.614** (gap −0.124 → −0.068),
  screening **0.521 → 0.625** (−0.254 → −0.149). K sweep is monotone to 3 and flat at 4. Nothing
  regresses beyond noise; the rare stratum mhcmatch already wins stays won. The gain is concentrated
  in **DP** (mean per-allele ΔAUPRC +0.108 vs DR +0.037) — DP scored 0.11–0.42 under a single PWM
  against DR's 0.6–0.94, so the human class-II "frequent gap" was largely a DP gap. See the benchmark
  repo's `bench/results/motif_mixture_mhc2.md`.
- **Calibrated MHC-II paths are ~3× slower** — this is where the mixture's cost lands, and only here.
  `restriction(calibrated=True)` per-peptide ~5.8s → ~17s; the `RankCalibrator` build ~17s → ~67s;
  `predict` likewise. The fast paths are untouched: default `restriction` (vote/enrichment, builds no
  `AnchorModel`) and `mhcmatch.ligand` span ranking (never calls `AnchorModel.score`). Set
  `n_motifs=1` to recover the previous speed. MHC-II model build 2.1s → ~19s (opt-in, once).

### Notes

- **What the components are not:** they come back 90–98% the *same* motif (per-anchor JS 0.02–0.05 of
  a possible 1.0), so this is not "each allele has two distinct binding motifs." Since `_m_step` gives
  each component its own best frame, the gain is plausibly a richer *register* model, not a richer
  motif model — recorded as untested. This also sidesteps the GibbsCluster multi-allele-deconvolution
  concern (its clusters are co-eluted *alleles*; our corpus is allele-labelled).
- **Measured on human MHC-II only.** Mouse and the interaction with the `%rank`/calibration accuracy
  are unvalidated; changing `n_motifs` back to 1 is the escape hatch.
- Doc fix: `load_markov1`'s docstring claimed `background="markov"` lifts MHC-I rare screening AUPRC
  ~+0.02; the committed tables say −0.019 (a sign flip). Corrected.

## [0.6.1] — 2026-07-17

### Fixed: the Potts affinity model's 8-mer encoding collision (code; weights deferred)

`PottsAffinity` encoded an MHC-I peptide as `core[:5] + core[-4:]`. For an 8-mer that puts index 4 in
two slots (`+5` and `−4` both land there), so the residue contributed two perfectly-correlated field
terms and a double-weighted coupling — the same defect v0.5.0 fixed for `AnchorModel` and never
propagated to the affinity head. Both the scorer and the trainer (`train_potts.py`) now route MHC-I
through `store.mhc1_positions`, the de-duplicated mapping. The two encodings agree for every L ≥ 9, so
only 8-mers were affected.

**The shipped weights are unchanged and 8-mer scores are unchanged** (bit-exact no-op, verified over
400 random 8–11mers). The encoding is bound to the weights by a version field in the `.npz` meta:
`PottsAffinity` uses the legacy slice for the shipped v0.6.0 weights and switches to the de-duplicated
mapping only for weights refit with it, so training and inference can never disagree about an 8-mer.
The numeric refit is deferred — the shipped `.npz` cannot currently be reproduced from the (gitignored,
regenerable) training data even with the legacy encoding, so a fresh fit would change every MHC-I
score for reasons unrelated to this defect. Tracked in the benchmark repo's
`results/potts_mhc1_encoding_defects.md`, which also documents the still-open length-blindness (defect 1).

## [0.6.0] — 2026-07-17

**MHC-II scoring changes by default**, and two gates that were measuring the wrong thing are fixed.
No API breaks; `AnchorModel(register="max")` restores the previous score.

- **MHC-II `score` integrates the binding register out** instead of maxing over frames. Every stratum
  × metric improves against NetMHCIIpan-4.3i; the rare stratum flips to winning all three. Frequent
  AUPRC gap −0.174 → −0.124.
- **The binder gate was a length detector** — a random 21-mer passed 98% of the time. Now a
  length-conditional `%rank ≤ 2`, MHC-II only; `restriction(cls="mhc1")` is byte-identical.
- **`predict_windows` was ~20× slower than it needed to be** — `_windows()` rebuilt an `AnchorModel`
  per binder (~10s each, ~20h over a 7,460-binder cohort) and re-derived the register from the wrong
  model, so the synthesised peptide could be cut from a frame the reported anchors did not describe.
- **The bench harness served stale examples** from a cache keyed on CLI args while the eligible
  allele set changed underneath. Caching is gone.
- **`bench/` now lives in [2026-mhcmatch-benchmark](https://github.com/antigenomics/2026-mhcmatch-benchmark)**;
  `bench/results/*.md` referenced below resolve there.

### Fixed

- **The MHC-II binder gate was a length detector.** `Store.restriction(diffuse=True)` gated on
  `anchor_score > 0.0`, but `AnchorModel.score` is a max over the `L−8` register frames, so it climbs
  with peptide length on **pure noise**: a random 15-mer was called a binder 85% of the time, a random
  21-mer **98%**. The gate now uses `percent_rank(allele, score, length=len(peptide)) <= 2` — a null
  of random peptides at the *query's own length*, so it goes through the same frame-max and the bias
  cancels (no independence assumption, unlike an extreme-value correction; overlapping frames are
  correlated). False-positive rate is now flat in length (3.7–6.7% for L=9…21) and is an explicit
  dial: `%rank <= t` passes `t%` of the null by construction. **Class-gated to MHC-II**: MHC-I is
  end-anchored with no frame max, and its length preference is real modelled biology that a
  length-conditional null would delete — `restriction(cls="mhc1")` is byte-identical and pays no
  calibration cost. Sensitivity on real held-out ligands goes 98% → 45% end-to-end; the old 98% was
  meaningless next to a 95% false-positive rate. No benchmark moves (`run_compare` scores
  `AnchorModel.score`, never `restriction`). See `bench/results/binder_gate_length_bias.md`.
- **The benchmark harness cached stale results.** `run_compare.py` keyed its `(examples, NetMHC
  scores)` pickle on the CLI args only — but `examples` depends on the eval-allele set, and
  `select_eval_alleles` gates on `a in pseudo`, so v0.5.0's pseudosequence fix silently changed which
  alleles are eligible while the key did not. The harness then served examples built from a **stale
  eval set** (rare n=21 against the true 24), producing numbers that disagreed with the committed
  results. **All disk caching is removed** from `run_compare.py`, `sample_concordance.py` and
  `bench/affinity/eval.py`; every run regenerates (a 35–70 s NetMHC sweep). The uncached harness now
  reproduces `compare_mhc1_human_hard_ligandbg.md` byte-identically.
### MHC-II scores now integrate the binding register out instead of maximising over it

`AnchorModel.score` for MHC-II was `max_r s_r` over every 9-mer core frame, which throws away *where*
the core sits. It now defaults to a marginal likelihood, `log Σ_r P(r | L, allele)·exp(s_r)`, under a
learned per-allele core-offset prior.

The prior is real signal, not bookkeeping. Real class-II cores sit ~3 residues from the N-terminus
(the groove protects the core while exopeptidases erode the flanks), so their offset distribution is
sharply peaked — DRB1_0101 15mers, H/Hmax **0.670** — while the *same model* lands uniformly on random
peptides (**0.998**). A decoy's argmax frame therefore sits at a low-prior offset about as often as
not while a real ligand's sits at the peak, and because the prior is normalized *within* a length the
term survives length-matched decoys rather than cancelling.

**Measured, head-to-head vs NetMHCIIpan-4.3i (seed 0, shortlist, identical examples): every stratum ×
metric improves and none regresses.**

| task | stratum | metric | `max` (old) | `marginal` (new) | Δ |
|---|---|---|---|---|---|
| allele-specificity | rare | AUPRC | 0.454 | **0.515** | +0.061 |
| allele-specificity | frequent | AUROC | 0.880 | **0.893** | +0.013 |
| allele-specificity | frequent | AUPRC | 0.508 | **0.557** | +0.049 |
| screening | rare | AUPRC | 0.555 | **0.652** | +0.097 |
| screening | rare | PPV@P | 0.376 | **0.541** | +0.165 |
| screening | frequent | AUPRC | 0.467 | **0.524** | +0.057 |

The rare stratum flips from losing AUPRC/PPV@P to winning all three metrics on both decoy modes (not
significant at n=19). The frequent AUPRC gap to NetMHCIIpan closes -0.174→-0.125 (hard) and
-0.308→-0.250 (screening) — narrowed, not closed.

Cross-allele ranking (`cv_mhc2_human_full.md`, 5-fold CV) improves too — top5 0.327 → **0.422**,
frequent recovery@5 0.298 → **0.409**, non-binder AUROC 0.556 → 0.596 — with **one exception**: rare
recovery@5 is flat-to-slightly-down (raw 0.490 → 0.487, diffuse 0.455 → 0.438), both inside one SD.
A rare allele has too few ligands to estimate its own offset shape, so it borrows one from groove
neighbours and there is little allele-specific offset signal left to add. Cross-allele diffusion
remains neutral-to-negative for MHC-II; this work does not change that.

- **Changed (MHC-II only):** `AnchorModel(register="marginal")` / `Store.anchor_model(register=...)`
  is the new default. Pass `register="max"` for the previous behaviour. MHC-I is untouched (it is
  end-anchored, so there is no register to integrate).
- **Unchanged:** `AnchorModel.best_register` still returns the argmax frame, so `decompose`, logos and
  the Potts affinity register oracle are unaffected. MBP85-99 / DRB1\*15:01 still ranks 2/149.
- **Cost:** MHC-II scoring 105k → **92k peptide-allele/s** (−12%; the prior is a cached per-(allele,
  length) lookup plus a logsumexp over frames that were computed anyway). Model fit is unchanged
  within noise (2.85s vs 2.86s on the 72k-peptide human shortlist panel) — the prior is estimated
  from the register-EM's existing frame assignments rather than a separate pass over the data.
- **Re-baselined:** `bench/results/register_em_mhc2.md`, `compare_mhc2_human_hard_ligandbg.md`,
  `compare_mhc2_human_random_proteomebg.md` — each keeps the old column alongside the new.
- **Does not fix the binder gate.** Marginalizing halves the length inflation (random peptides,
  9mer → 21mer: +4.44 nats → **+2.28**) but leaves a Jensen residual, so a random 21-mer would still
  pass a raw-score gate two thirds of the time. The gate is fixed separately and orthogonally by the
  length-conditional `%rank` above.

### Assay provenance: the panel is not what SOURCES said, and the benchmark can now say so

`bench/affinity/SOURCES.md` claimed the presentation tables "keep eluted-ligand positives only".
**False** — **36,881** class-II (epitope, allele) pairs have no mass-spectrometry assay at all
(14,969 competitive-radioactivity, 13,416 high-throughput multiplexed, 8,343
competitive-fluorescence, 237 Edman degradation). What the tables drop is the quantitative
*measurement*, not the binding-assay *rows*. Both SOURCES files are corrected.

New: `bench/compare/provenance.py` + `run_compare.py --el-only`, an **evaluation stratum** that makes
only mass-spec-supported pairs eligible as positives. **Training still uses the whole corpus** —
binding-assay peptides do bind, so they are valid motif evidence, and the house rule is one corpus
tuned per task by parameter (`CLAUDE.md`), never a smaller training set to make a benchmark look
clean. Assay type is absent from the pmhc schema, so it is joined from the raw IEDB dump on
`(epitope, reference_id)` — present in both tables, so no restriction-name parsing — and cached
(3.19M pairs, ~90s to build).

**Source-conditioning was tested and rejected.** The obvious refinement is an adjusted general model
per provenance, since EL boundaries are biological (offset H/Hmax 0.720) and binding-assay boundaries
are experimenter-chosen (0.990, flat as random). Held out, the corpus-learned offset prior beats a
uniform one by **+0.010** on EL queries and **+0.001** on BA queries — it helps where boundaries carry
information and is harmless where they do not. The general model already serves EL, BA and in-silico
queries; no `source` switch is warranted.

**The share is confounded with allele, which is what makes it matter:**

| panel | frequent alleles | thin alleles | alleles with zero EL |
|---|---|---|---|
| human class II | 25.7% non-MS | 83.1% non-MS | **15 of 52** |
| mouse class II | H-2-IAb 4% non-MS | H-2-IEd/IAs/IAq ~100% | **6 of 13** |

- **The human `rare` stratum has no eluted-ligand positives to evaluate on** — 15 of 52 alleles have
  zero eluted ligands, 8 more are under a 20-ligand floor. mhcmatch's rare-stratum win
  (`compare_mhc2_human_hard_ligandbg.md`, AUROC 0.842 vs 0.813) therefore answers "reproduce IEDB",
  not "find eluted ligands". Both are real questions; the pair is reported.
- **It does not move the gap.** Both tools score higher on eluted-ligand positives, and the frequent
  gap barely shifts (AUROC -0.053 → -0.050, AUPRC -0.124 → -0.124). It changes what a number is
  *about*, not who wins.
- Binding-assay rows stay in training — those peptides do bind, so they are valid *motif* evidence.
  What they are not is evidence about *boundaries* (`bench/results/length_prior_mhc2.md`).

### First mouse MHC-II head-to-head — two tables, two questions

Both are reported; neither supersedes the other.

**`compare_mhc2_mouse_hard_ligandbg.md` — reproduce IEDB's mouse annotation. mhcmatch wins all nine
cells**, the only panel where it leads every stratum on every metric (medium AUROC +0.422,
AUPRC +0.424, p<0.001). Recorded observation: NetMHCIIpan's medium AUROC is 0.464, below chance —
mouse provenance is confounded with allele (H-2-IAb 96% mass-spec over 10,797 peptides; H-2-IEd/IAs/IAq
0%), so a BA-only allele's positives face I-Ab's real-ligand decoys and an EL-trained tool ranks the
decoys higher. `n` is 1/4/3 alleles of 13.

**`compare_mhc2_mouse_random_proteomebg.md` — find eluted ligands (`--el-only`, proteome decoys).**
NetMHCIIpan is above chance everywhere and nothing separates the tools: AUROC 0.793 vs 0.789
(+0.004, p=0.94), NetMHCIIpan's AUPRC lead inside its own interval (0.256 vs 0.320, p=0.49). Three
alleles — H-2-IAb (7,990 EL), H-2-IAd (161), H-2-IEk (97) — of a 13-allele panel.

This does refute the idea that mouse is the "uncontaminated axis" — the obstacle was never
NetMHCIIpan's thin mouse training, it is the panel's provenance imbalance.

- **Fixed:** `run_compare.py` hardcoded `human.fasta.gz` as the decoy proteome regardless of
  `--species`. Measured impact was small (KL(mouse‖human) over proteome AA frequencies = 0.00043
  nats), but the flag was being ignored. `PROTEOME_AA_FREQ` / `proteome_markov1.tsv` stay human as a
  documented approximation.
- `provenance.el_only(min_peptides=20)` drops alleles too thin to support a metric, and **logs** what
  it dropped. Without the floor the mouse "rare" stratum is three alleles with 2, 3 and 11 ligands,
  where mhcmatch "wins" AUROC by +0.248 and the opponent's PPV@P is a coin flip.
- **`predict_windows` synthesised the wrong register (MHC-II).** `_windows()` called
  `store.anchor_model("mhc2")` with the *defaults* (`footprint="anchor"`, `background="ligand"`) — a
  different model from the `adaptive`/`proteome` one that had just scored the peptide — and re-derived
  the binding register from it. So `synth_peptide` / `model_peptide` could be cut from a different
  register than the one `anchors` / `tcr_facing` / `agretopicity` were reported for, breaking the
  invariant asserted in the comment directly above the call. The scored register was already in scope
  and is now passed in. `synth_peptide` is what gets ordered as a peptide, so this was a correctness
  bug, not a cosmetic one.
- **The same call rebuilt an `AnchorModel` per binder.** An MHC-II `AnchorModel` costs ~10 s to build
  and `_windows()` ran once per kept binder — ~20 h of rebuilds over a 7,460-binder cohort. Passing
  the register in removes the call entirely.
- **`build_scorer` is now memoised on the store.** It depends only on the panel, never on the query
  alleles, so scoring many samples against one store reuses a single build instead of paying the
  MHC-II model and calibrator per call. Measured on a real sample: 39.6 s cold → 0.0 s warm.
- **`agretopicity` was computed from the rounded WT nM.** It divided the unrounded mutant IC50 by
  `wt_affinity_nm`, which is rounded to 1dp for display, while `dai` recomputes both unrounded — so
  the two disagreed by up to ~0.5% and could report opposite directions for the same peptide near
  agretopicity 1. Now divides the unrounded pair (the displayed field keeps its rounding). The
  `amplitude` field comment also claimed `A = Kd_WT/Kd_MT`, omitting the saturation correction
  `affinity.py` applies — which reads as "amplitude == 1/agretopicity", and it is not.
- **`bench/compare/sample_concordance.py` read the class-II pipeline column with the sign flipped.**
  The pipeline renames TLimmuno2's `Rank` to `affinity`, so it is a rank fraction (lower = stronger,
  gated at 0.1), not TLimmuno2's `prediction` (higher = stronger). It negates like class I.
  `score_epitopes.py` had it right; the bench reader did not. Part of why
  `bench/results/concordance_tesla1_mhc2.md` reports mhcmatch~pipeline = −0.034.

## [0.5.0] — 2026-07-16

**Allele coverage was broken: 68% of MHC-I and 80% of MHC-II alleles could not be resolved at all.**
Plus the MHC-I score becomes length-aware by default. No API breaks; some defaults change (below).

### Fixed

- **Pseudosequence name index (the headline).** Alleles sharing a 34-mer groove collapse to one FASTA
  record, but only the *first* allele's name was written — the other **8854 of MHC-I's 12997** and
  **8839 of MHC-II's 11048** were silently unresolvable. Not rare variants: `HLA-B*14:02`, `B*18:05`,
  `C*03:04`, `C*03:02` all returned nothing while `HLA-C03:438` shipped. `restriction()` and
  `predict()` gave no answer for any of them. The collapse was always right; the name index was lost.
  Headers now list every allele of the group; each resolves to **its own true 34-mer** (the group is
  exact-identity, so this is not a nearest-neighbour guess).
- **MHC-I 8-mer anchor collision.** `MHC1_CORE`'s `+5` and `−4` both mapped to index 4 of an 8-mer,
  double-counting it in the score *and* filing one residue under two positions during training.
  `store.mhc1_positions` is now the single de-duplicated mapping shared by scorer and estimator.
  **8-mer scores change.**

### Added

- **IPD-IMGT/HLA as a second pseudosequence source** — **+7085 class-I alleles** (20082 total, 5407
  unique grooves). NetMHCpan's table lags IMGT and omits **HLA-F entirely**. The 34 positions are
  recovered from the alleles the table already covers, cross-checked between genes (HLA-B and HLA-C
  solve independently and agree), and verified by re-deriving every known allele: **21935 exact, 4
  mismatch (0.018%)**. NetMHCpan wins every conflict, so no covered allele changes. The human MHC-I
  reference panel goes **166/203 → 203/203** scorable. Regenerate with `bench/build_pseudo_fasta.py`
  (now vendored here; mhcmatch no longer re-syncs this data from `tcren`).
- **DP/DQ α-chain imputation for lookup** (`pseudoseq.alpha_prior`, `data/mhc2_alpha_prior.tsv`).
  MHC-II is an αβ heterodimer but 1.5% of panel records type only β. `HLA-DPB1*11:01` returned `nan`;
  it now resolves to `HLA-DPA10201-DPB11101`. Learned from the panel, keyed on **P(34-mer groove | β)
  ≥ 0.95 over ≥ 50 ligands** — the groove, not the allele name or its 2-digit group (`DQA1*01:02` and
  `DQA1*01:05` share the group but not the 34-mer). Rediscovers DQ2.5 and DQ8 from linkage
  disequilibrium. 9 rare DQ βs fail the bar and stay unresolved on purpose.

### Changed

- **`length_prior` and `length_motifs` now default ON for MHC-I.** The anchor log-odds summed a
  length-invariant number of terms, so a 10-mer and a 9-mer with the same anchors scored
  bit-identically — while a length-only classifier reaches maxF1 0.802 on the MixMHCpred3 benchmark.
  Adds a per-allele ligand-length factor (kernel-shrunk over groove pseudosequences, so rare alleles
  borrow a length profile from neighbours) plus per-length motifs with an exact backoff: an allele
  with no ligands at length L reproduces the pooled model bit-for-bit and provably cannot regress.
  MHC-II is untouched (both are class-gated). Pass `length_prior=False, length_motifs=False` for the
  old behaviour. Costs ~9% throughput.
- **`Store.from_records`/`from_pmhc` gain `impute_alpha` (default OFF).** Opposite to the lookup path,
  and measured: admitting β-only records to the reference *panel* moves held-out AUROC −0.0019 and
  AUPRC −0.0012 over the 13 affected alleles, worst where the merge is biggest (`DPB1*11:01` +89%
  ligands → −0.0155 AUROC). A study that skipped α-typing produced noisier ligand calls too.

### Benchmarks

MixMHCpred3 (20 HLA-typed samples, leak-free panel; MixMHCpred3.0 = 0.911, BigMHC = 0.911,
NetMHCpan4.1 = 0.899):

| | maxF1 |
|---|---|
| 0.4.2 | 0.8501 |
| **0.5.0** | **0.8907** |

Length work +0.0306 and the name-index fix +0.0104 are additive (+0.0410 predicted, +0.0407 measured).
The IMGT source is worth **0.000 here by design** — every benchmark allele was already covered; it buys
coverage, not score. `bench/results/compare_*.md` are regenerated.

**The head-to-head numbers moved and the eval set moved with them** — `select_eval_alleles` gates on
`a in pseudo`, so fixing the name index made previously-invisible alleles eligible (MHC-I rare 21 → 24,
MHC-II 37 → 47 total). The strata are **not comparable to 0.4.2's**, and NetMHCpan/NetMHCIIpan — fixed
binaries — moved too (MHC-I rare AUROC 0.971 → 0.945; MHC-II rare 0.858 → 0.881), which only the eval
set changing can explain.

- **MHC-I allele-specificity improved**: rare went from −0.021 AUROC (NetMHCpan's) to **+0.008** (a
  wash); frequent AUPRC 0.812 → **0.850**. Medium/frequent stay significant wins (p < 0.001).
- **MHC-II**: on a *frozen* eval set the model change alone is +0.0008 AUROC / −0.0107 AUPRC — and that
  AUPRC delta is one allele with a **single** ligand (`DRB1_0302`, held out, hence scored zero-shot)
  moving one rank. 95% CI [−0.0367, +0.0029], 31/40 alleles same-or-better, frequent stratum +0.0002.
  No regression.

## [0.3.0] — 2026-07-14

**Core → full presented ligand** (`mhcmatch.ligand`), plus the register refactor it needed. Backward
compatible: new module and one new `AnchorModel` method; existing defaults unchanged.

### Added

- **`mhcmatch.ligand`** — extend a 9-mer MHC-II binding core to the peptide that is actually presented.
  Three evidence tiers: `observed` (a real eluted ligand containing the core), `modeled`
  (`SpanModel`, a flank/context model fit to mass-spec ligandome data), `fixed` (caller flanks,
  clipped at protein bounds and *reported*, never silently shortened).
  - **Not a cleavage predictor, by design.** MHC-II is bind-first-trim-later, so there is no strong
    sequence-specific endoprotease step to simulate; the one dedicated MHC-II cleavage motif
    (PMID 30127785) gets AUC 0.767 on ligands and *zero* on CD4 epitopes. The model is the learned
    flank model the field actually uses (NetMHCIIpan `-context`, PMID 30446001; MHCflurry-2.0
    processing, PMID 32711842): 12 terminus-relative context positions vs an order-1 Markov proteome
    null, plus a ligand-length prior. Allele-agnostic (measured: per-allele JSD 0.003–0.010), **no
    free parameters**.
  - **Not an immunogenicity predictor** — context is documented to *degrade* CD4 epitope benchmarks
    (PMID 32406916). It answers "what ligand?", not "is it immunogenic?".
  - `processing_score()` for MHC-I (the peptide *is* the ligand, so it returns a score, never a span).
    Class I and class II are separate entry points with **no class inference** — a 9-mer class-II core
    is always ≤11 and would misroute.
  - **`STRUCTURE_FLANK = 2`** (13mer) and **`ASSAY_FLANK = 6`** (21mer) — the recommended fixed flanks,
    both measured. The span model's point estimate is *not* accurate enough to pick a peptide from
    (both boundaries within ±2 only 47% of the time, barely beating a centred 15mer), so these are the
    defaults to use; the model answers "what was eluted?", not "what should I make?".
- **`AnchorModel.best_register(peptide, allele) -> (start, score)`** — returns the winning register
  frame that `score()` already computed and discarded. `score()` and `_refit_registers()` now collapse
  onto it (bit-identical). The three heuristic-register duplicates collapse onto `store._mhc2_register`.
  The two registers stay two **by design** (ROADMAP §7).
- **`mhcmatch span`** CLI subcommand.
- `bench/train_spans.py`, `bench/bench_spans.py`, `bench/pdb_flanks.py`;
  `bench/results/spans_mhc{1,2}_human.md`; vendored `data/ligand_context.tsv`.

### Fixed / found

- **Documented an open bug: the MHC-II binder gate is a length detector.**
  `Store.restriction(diffuse=True)` gates on `anchor_score > 0.0`, but `AnchorModel.score` is a max
  over register frames and grows with length even on noise — a **random 21-mer passes 98%** of the
  time, a random 15-mer 85%. Not fixed here (it changes `restriction()` semantics);
  `bench/results/binder_gate_length_bias.md`.

### Measured

- MHC-II span recovery (gene-split, leak-canaried): set-recall **0.158** vs 0.069 for centring a
  15mer, against a **0.547** nested-set oracle ceiling. Honest caveat: it does *not* beat that
  baseline on mean boundary error.
- MHC-I context: full 12-position AUROC 0.814, but **flank-only (the honest processing signal) 0.558**,
  shuffled control 0.501.
- 93 real pMHC-II crystals (Canonical2026): resolved peptide **median length 13**, median 2 ordered
  flanking residues per side; only **13%** resolve ≤11 residues — so core±1 is too short.
- Known-biology control: Pro **2.00×** enriched inside the ligand, **0.25×** depleted in the flank
  (the aminopeptidase stop signal) — the *opposite* sign to the naive prior.

## [0.2.0] — 2026-07-09

First head-to-head against NetMHCpan, plus the scoring and reporting upgrades it motivated. All
additions are backward-compatible (new opt-in parameters; existing defaults unchanged).

### Added

- **Head-to-head benchmark harness** (`bench/compare/`) vs **NetMHCpan-4.2b** / **NetMHCIIpan-4.3i**
  on two shared per-(peptide, allele) tasks — *allele-specificity* (decoys = other alleles' ligands)
  and *presented-vs-random screening* — stratified rare/medium/frequent, with AUROC / **AUPRC / PPV@k**,
  bootstrap CIs and paired **DeLong**/bootstrap significance. Results in `bench/results/compare_*.md`,
  provenance in `bench/compare/SOURCES.md`. Caches `(examples, NetMHC scores)` for fast model iteration.
- **Calibrated outputs** (`mhcmatch.calibrate`): per-allele **%rank** vs a random-peptide background
  (NetMHCpan `%Rank_EL` analogue), isotonic **P(present)**, and a qualitative binding **band**. Wired
  into `Store.restriction(calibrated=True)` and the CLI (`mhcmatch restriction --calibrated`).
- **`AnchorModel` scoring footprints** (`footprint=`): `"anchor"` (default, primary pockets),
  `"core"` (full binding core), `"adaptive"` (class-aware — anchors for rare MHC-I alleles, full core
  for MHC-II and well-sampled MHC-I).
- **`AnchorModel` log-odds nulls** (`background=`): `"ligand"` (default, allele-*specificity*),
  `"proteome"` (presentation — `log(θ_A / p_proteome)`, far better at ligand-vs-random screening),
  `"markov"` (order-1 proteome conditional, a rare-allele lift). New vendored
  `data/proteome_markov1.tsv`.

### Results (shortlist tier, human, seed 0)

- **Allele-specificity:** mhcmatch **beats** NetMHCpan on MHC-I medium+frequent (AUROC/AUPRC/PPV@k,
  p<0.001; frequent AUPRC 0.81 vs 0.69).
- **Screening (proteome null):** mhcmatch **beats** NetMHCpan on MHC-I medium+frequent AUPRC/AUROC and
  NetMHCIIpan on MHC-II rare AUPRC (0.69 vs 0.58). Rare MHC-I remains NetMHCpan's.
- **Speed:** ~68× faster than NetMHCpan (195k vs 2.9k peptide-allele scores/s).

## [0.1.0]

Initial release — restriction/presentation, similarity search, anchor/TCR-facing split, source
lookup, motif logos, and the pseudosequence cross-allele diffusion model.
