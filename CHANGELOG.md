# Changelog

All notable changes to `mhcmatch`. Format loosely follows [Keep a Changelog](https://keepachangelog.com);
versioning is [SemVer](https://semver.org).

> Note: 0.4.0–0.4.2 shipped without entries here. This file jumps 0.3.0 → 0.5.0; see `git log` for
> the 0.4.x range.

## [1.0.4] --- 2026-08-24

The command line can now emit every table a figure needs.

### Added

- **`mhcmatch rank --coefficients` and `mhcmatch rank --holdout`.** The fitted EPIC aggregate had no
  CLI equivalent, so anyone wanting the model itself --- rather than a scored candidate list --- had
  to import the package and read `data/aggregate_mhc1.json` by hand. `--coefficients` prints it as
  TSV (block, term, coefficient, Laplace and bootstrap sd, *z*, *p*, the 95 % cluster-bootstrap
  interval, sign stability); `--holdout` prints the same fit's nine leave-one-screen-out AUROCs with
  their decided/undecided flag, plus both grouped cross-validations. Both read the shipped artifact
  and refit nothing, so a figure built on them and a run of `rank` are the same model by
  construction rather than by a comparison someone has to remember to make. `mode` and `input` are
  optional under either flag.
- **`mhcmatch rank pairs FILE`** --- a third input shape, a TSV of `peptide` / `wt_peptide` /
  `allele` (+ optional `gene`, `tpm`). `rank fasta` needs mutation-spanning windows and `rank table`
  needs another tool's `.scored.csv`; neither is what you have when a caller has already given you
  the mutant k-mer, its germline counterpart and the restricting allele --- which is how every
  neoantigen screen is distributed. Scoring one therefore meant reimplementing `rank` outside the
  package, and a reimplementation is a second model nobody benchmarked. Rows are grouped by allele
  and each group scored in one `predict.binder_ranks` call, so the per-allele calibrator background
  is paid once per allele rather than once per row; the wild types go through the same call, so the
  ratio between a WT and a mutant IC50 is a property of the substitution and not of two code paths.
  A row with no wild type is kept: `wt_absent` carries it and agretopicity stays undefined rather
  than zero. Exposed as `mhcmatch.rank.rank_pairs`.
- **TSV output for `scan`, `logo` and `expression`,** under `--out FILE` or `--tsv`. All three
  printed aligned text only. `expression` wrote `median 0.33` and `IQR 0.1-0.9` *inside* cells,
  which reads well and parses badly; `logo` kept the top three residues per position where the TSV
  now carries the whole PWM; `scan` collapsed a window's alleles into one comma-joined cell, and the
  TSV gives one row per (window, allele) with its enrichment and vote count. The aligned form is
  unchanged and is still the default.

### Changed

- **`cassette score --pool` could not read the pool `cassette select` had just written.** It
  matched each unit to a pool row by comparing scores to 1e-9, and `select` writes six decimal
  places where a pool written at six significant figures does not survive an exact float
  comparison --- so the one chain `docs/cassette.rst` recommends reported the cassette's own units
  as absent from the pool they came from. Where both sides carry peptides the peptide is now the
  key, which is what a unit's identity actually is; the score stays the fallback for a caller
  passing bare vectors.
- `expression --safety` no longer requires `--tissue`/`--tumor`. `lookup()` demands exactly one of
  them, and `cmd_expression` called it before reaching the safety block, so asking only for the
  cross-tissue profile raised instead of answering.
- Shipped artifacts re-stamped to 1.0.4 (`mhcmatch build`): three anchor models and
  `corpus_tables.npz`. Contents unchanged; `mhcmatch build --check` reports 0 stale of 27.

## [Unreleased]

Documentation only.

### Added

- **The charge column is named where it is used.** `C_phys_charge` is half of EPIC v4's chemistry
  block and the README described that block as if burial were alone in it --- it said
  `C_phys_buried` was "the chemistry block's only fitted term", which the shipped artifact
  contradicts. The README now carries the arm: Atchley 2005's fifth factor, electrostatic charge,
  averaged over the TCR face, selected on its **residual against Rose** rather than on its own
  AUROC, orthogonal to burial at *r* = **+0.008** per peptide against **−0.837** for the Kidera KF4
  pair it replaced, and paying for itself in burial's stability --- bootstrap sd **0.0874 →
  0.0487**, *z* **+1.71 → +2.34**, *p* **0.088 → 0.020**, sign stability **96.5 % → 100 %** over 400
  cluster bootstraps on 354,909 rows and 958 immunogenic peptides. Its own coefficient stays what it
  is, **−0.0634 (z −1.21, p = 0.225)**: the column being fixed is burial.
- **What MJ1996 and TCRen are for, in `docs/complementarity.rst` and the README.** They are 4 of the
  30 columns of `complement.score` and are in **no** EPIC v4 term, so what they carry is a physics
  distinction rather than a ranking. Stated with the measurement behind it: MJ1996 is **96.4 %
  one-body** with its two leading modes at exactly **±0.851** against Kyte–Doolittle, TCRen is
  **3.29 %** --- below its own composition-matched shuffle floor of **9.68 ± 2.16 %** over 500
  shuffles --- with neither leading mode a hydropathy axis. The face assignment is checked on
  coordinates: on 3,875 (structure, position) rows from 374 crystals, MJ1996 alone detects
  groove-floor contact at AUROC **0.5818** and the TCRen marginal at **0.4801**, below chance, which
  is the ordering the block assumes; class II reproduces neither (**0.5171** / **0.5368**, 94
  structures) and the docs say so. Fitted coefficients per column standard deviation are recorded
  alongside --- `mj_anchor` **+0.0354**, `mj_tcr` **−0.0927**, `para_tcr` **−0.0299**,
  `para_sd_tcr` **+0.0052**, against `aa_tcr`'s **−0.8796**.
- **[`antigenomics/tcren`](https://github.com/antigenomics/tcren) is credited and linked** from the
  README, `docs/complementarity.rst` and `docs/property_basis.rst` --- the potentials, the
  374-crystal contact maps and the spectral analysis are upstream work there. The tables are
  vendored, so `tcren` stays a runtime dependency of the optional `[structure]` extra alone.

### Fixed

- **AF5's cysteine loading was transcribed wrongly in `docs/burial.rst`**, as **+0.0056**. That is
  `PROTFP:ProtFP8`'s value --- the row above it in `bench/results/physchem_pc.md`. AF5's is
  **−0.0028**, the lowest of the 141 complete scales swept, which is what
  `complement.PHYS_SCALE_CHARGE`, `rank.PHYS_COLUMNS` and the theory appendix already carried.
- **`skills/mhcmatch/SKILL.md` said `burial` sums over the TCR face.** It averages ---
  `per_residue=True` is the default and the same file says so 86 lines later, because the summed
  form was Pearson **+0.954** with peptide length. The cell now also names `C_phys_charge`.

## [1.0.3] - 2026-08-24

Documentation only. **No executable code changed** except the version string: six of the seven
touched modules are AST-identical to 1.0.2 with docstrings excluded, and the seventh is the bump.

### Changed

- **The user-facing docs stop being a changelog.** README, the Sphinx pages, the agent skill and the
  module docstrings carried a running account of models that no longer ship --- `BOECRT`, `PADEC`,
  `BDEVF`, `GRAND`, `ipred` --- and of when each behaviour arrived ("Since 0.20.0...", "removed in
  0.22.0", "the pre-0.19.0 ordering"). A reader installing the package now needs to know what it
  does, not what it used to do. Every one of those is gone from the README (44 KB -> 39 KB), the
  docs, `skills/mhcmatch/SKILL.md` and the notebooks; the history lives here and in the benchmark
  repository's `MODELS.md`, which is where it belongs.
- **`docs/models.rst` removed.** It was a naming registry for retired models. Its one live section,
  the occupancy-versus-agretopicity derivation, moved to `docs/neoantigen.rst` and keeps its
  `occupancy-vs-agretopicity` label so existing cross-references resolve.
- `docs/complementarity.rst` loses its 463-line `ipred` archive (1,109 -> 644 lines).

### Notes

- 1.0.1 and 1.0.2 were tagged and neither reached PyPI; the 1.0.2 publish run was cancelled during
  its artifact-regeneration step, before anything was uploaded. 1.0.3 is the first 1.x on PyPI.

## [1.0.2] - 2026-08-24

A cleanup pass over 1.0.1, which was tagged but never published. 1.0.2 was likewise
tagged and never published --- see 1.0.3. No behaviour changes: every
refactor below was checked bit-identical against the code it replaces before it was kept.

### Changed

- **One bisection, in one place.** `rank.probability` re-implemented `cassette.prob_offset`'s
  root-find character for character; it now calls it. Verified bit-identical on 10, 500 and 5,000
  scores including non-finite entries. The difference between the two functions is which batch you
  hand it, not the arithmetic, and two copies of the arithmetic is how that distinction gets blurred.
- **One clipped sigmoid in `cassette`**, `_p(scores, offset)`, replacing five hand-written copies of
  `1/(1+exp(-clip(s+b, -60, 60)))`. `offset` broadcasts, which is what `group_offsets` needs.
  `group_offsets` verified bit-identical on 3, 7 and 46 groups.
- **One text opener in the CLI**, `_open_text`, replacing four copies of the gzip-or-plain /
  stdin-or-file / `try`-`finally`-that-must-not-close-stdin dance. Parsing and error messages stay
  in each reader, because those are the part that differs and the part a caller reads when their
  table is wrong. Net -2 lines; the point is that the "never close stdin" invariant now has one home.

### Added

- `tests/test_structure.py`. `StructureScorer.__init__` calls `_require_tcren` before it does
  anything, so no scorer behaviour is reachable in a default install and this module's coverage is
  low **by construction**. The tests record which half is which: the vendored template table is
  well-formed and actually read, the missing-extra path names `mhcmatch[structure]`, and the scorer
  itself is `skipif`-gated on the extra.

### Packaging

- `Development Status :: 5 - Production/Stable` (it said Beta on a 1.0 release), explicit 3.10-3.13
  and `OS Independent` classifiers, and a real `[project.urls]` -- Homepage, Documentation,
  Repository, Changelog, Issues. There were two `[project.urls]` tables; the second silently won.
- `viz = ["networkx"]` keeps its place but says what it is for: **the library does not import
  networkx**. It is what the benchmark repository's promiscuity-graph figure generator needs.

### Considered and not done

- **`mhcmatch.recognition` is not dead code.** It has no caller inside the package, which is what an
  import graph shows, but `bench/neoag/features_grand.py` calls `recognition.score` to build the
  **EPIC fitting frame** and the theory appendix describes it as shipped. Removing it would break the
  shipped model's own provenance.
- `mhcmatch.luksza` stays. `EPIC` retired the term in 0.21.0 and only the `__init__` re-export
  imports it, but it is the published Luksza quantity and computing it should not require the
  benchmark repository.
- The README is 44 KB and is the PyPI long description. That is a lot, and it is the author's prose;
  trimming it is an editorial decision, not a refactor.

## [1.0.1] - 2026-08-24

First 1.x. The version skips 1.0.0 deliberately (`ROADMAP.md` sec. 5e).
**Tagged but never published to PyPI** --- superseded by 1.0.2 before the release gate opened.

### Added

- **`mhcmatch.cassette` — the two operations a cassette needs, as one module.** Everything below
  existed as an algorithm somewhere in the benchmark repository; what it did not have was a shipped
  home, a CLI, or a test.

  | | |
  |---|---|
  | `select(scores, peptides, alleles, k, tol)` | choose *k* units (+/- `tol`) maximising `H = sum h - sum J`, the mean-variance objective derived from the design goal. Greedy `O(kN)` plus a bounded swap pass; reaches the brute-force optimum on every pool small enough to enumerate |
  | `score(...)` | expected responding units, `P(>= k)` under the block model, `n_effective`, allotype coverage, the three pairwise statistics, and `lam` |
  | `lam(h, sel, k)` | **the axis that compares cassettes across donors and across sizes** — `H(S)` minus the exact log partition function over every size-*k* subset of that donor's own pool, plus `log C(N, k)`. Zero is a uniformly random subset of the same pool |
  | `prob_offset` / `group_offsets` | the calibration offset, fitted over a batch that does not move, or one per group. **These report different quantities** — a level and an enrichment — and the docstrings say which |
  | `goal_energy`, `greedy`, `refine`, `overlap`, `pair_stats`, `log_ek`, `energy` | the primitives, each with its closed form checked against the `O(k^2)` or brute-force sum it replaces |

- **`mhcmatch cassette` — the CLI, with sub-verbs.** `select`, `score`, `build`, `order`, `deslip`.
  The first two-level command in this CLI; the `-v`/`-q` loop in `main` descends into it, because
  otherwise `cassette select -v` is an unrecognised argument while `cassette -v select` works.

- **`MHCMATCH_CASSETTE_SCORE`, a Nextflow process that is deliberately not per sample.** It collects
  every donor's table and fits **one** calibration offset over the run. `rank` anchors `p_response`
  on the batch it is handed, so a per-donor call pins every donor's mean candidate probability to
  the declared prevalence: measured on 7,261 TCGA donors with pools spanning 1 to 5,221 candidates,
  every per-donor-anchored pool mean lands on **0.060163**, standard deviation **2.75e-17**. Two
  donors' numbers are then the same number. `params.mhcmatch_cassette_per_donor_offset` asks for the
  enrichment reading instead, and says in its own help what that costs.

- `portfolio.betabinom_rho(..., profile=False)` maximises over `(p, rho)` jointly and reports the
  fitted `p`. The benchmark repository had grown a second implementation of this estimator for
  exactly that; it now calls this one.

### Changed

- **`mhcmatch vector` is now `mhcmatch cassette build`**, and `mhcmatch deslip` is
  `mhcmatch cassette deslip`. Both old names remain as aliases for one release and print a
  deprecation line to stderr. `MHCMATCH_VECTOR` is `MHCMATCH_CASSETTE` in the Nextflow module; its
  `params.mhcmatch_vector_*` names are **unchanged**, because an unknown Nextflow parameter is
  ignored rather than rejected and a rename would silently drop every deployed config's settings.
- `cassette order` runs the assembly half alone on units already chosen, so `--n0` is not required
  there. Same code path as `build`, so the two cannot diverge.

### Notes

- `score` does **not** report `H`. `goal_energy` renormalises the overlap to the set it is handed
  and `overlap`'s dominance channel is scaled by that set's range, so an `H` computed on a cassette
  alone is not the `H` `select` maximised over the pool — and a rule that spent expected count on
  non-overlapping units would score identically to one that did not. To compare two rules on the
  objective, build `(h, J)` once over the pool and evaluate both index sets with `energy`. `lam`
  needs none of that and is the axis that already crosses donors and sizes.

## [0.27.0] - 2026-08-23

### Changed

- **The shipped scorer is EPIC v4.** `data/aggregate_mhc1.json` carries `"version": 4`, and four of
  the nine terms are respecified rather than refitted:

  | v3 | v4 | why |
  |---|---|---|
  | `binder` | `pres` | `binder` Fisher-combined the presentation `%rank` with the affinity rank, and `occupancy` already carries affinity at Spearman −1.000000 against `kd_mt`. `pres` is the presentation rank alone. |
  | `expr` + `expr_missing` | `expr_pct` | two terms to one: the within-cohort expression percentile, unit-free and needing no missingness flag (above) |
  | `C_phys_rose` | `C_phys_buried` | rename; the Rose 1985 scale *is* mean fractional area loss |
  | `C_phys_hydrop` | `C_phys_charge` | Kidera KF4 is burial measured a second way (r = −0.837 per peptide) and the pair was not identified. Atchley AF5 is orthogonal by measurement, r = +0.008. |
  | corpus kernel | Hamming → BLOSUM62 | identity-normalised, `K[u,x] = exp(κ(S[u,x] − S[u,u]))`, k = 3, sliced face |

  On 354,909 rows and 958 immunogenic peptides over nine screens: BIC 4212.2 → **4168.6**,
  leave-one-screen-out mean AUROC over the seven deciding screens 0.6503 → **0.6654** and median
  0.6325 → **0.6399**, twin-grouped five-fold cross-validation over the whole database 0.6325 →
  **0.6399**. Seven of the nine screens rank higher held out than under v3, by +0.0007
  (IEDB_neoag, 1,280 rows / 601 positives) to +0.0443 (VACCIMEL, 93 rows / 27 positives). ITSNdb
  falls 0.6566 → 0.6399 on 149 rows with 89 positives, and NCI, which has six held-out positives
  and does not decide, falls 0.8513 → 0.8248.

### Removed

- **EPIC v3 is gone from the library — every path, name and fallback.** It was carried as a
  parallel vocabulary so one library could score either artifact, and the cost of that was a
  second definition of every changed term, live, with nothing exercising it. What went:

  | removed | why |
  |---|---|
  | `C_phys_rose`, `C_phys_hydrop` from `AGGREGATE_COLUMNS` / `PHYS_COLUMNS` | computed and emitted on every row for a model nothing ships |
  | `binder` from the feature union `rank._finish` supplies | v3's presentation term; `pres` is the fitted one. `predict.binder_score` and the `binder` output column are unaffected — they are a prediction API, not a model version |
  | `complement.PHYS_SCALE_HYDROP` | replaced by `PHYS_SCALE_CHARGE = "ATCHLEY:AF5"`, the scale actually fitted |
  | the Hamming default in `mimicry.corpus_geometry` | an artifact naming no kernel was silently contracted against Hamming — a `kappa` fitted against a graded kernel scored that way is a different feature, not a smaller effect. It now raises. So does an unknown face mask; the `"tcr5"` alias is gone |
  | `supersedes` and `phys_scale_hydrop` from the artifact | v3 pointers, and `epic_v4_fit.py` no longer writes them |

  Four tests asserted the v3 half of a pair and now assert the v4 fact instead;
  `test_c_phys_buried_is_the_same_column_as_c_phys_rose` had nothing left to say and is deleted.

### Changed

- **The expression block is one term, `expr_pct`, and it is a rank.** It was `expr` (global *z* of
  `ln(1+TPM)`, mean-imputed) plus `expr_missing`, a binary indicator. The indicator was metadata:
  `expr_source` is very nearly a screen label — Neopep, 87% of the corpus, is 95.5% GTEx
  matched-tissue; NCI and HiTIDE are 100% measured; four small screens are 93–99% reference — so
  inside every screen but one the flag is a constant the per-screen intercept already carries.
  Measured by drop-one: ΔBIC +36.6, the second largest of the nine, against a held-out cost of
  0.6654 → 0.6624, the smallest of the four presentation and expression terms.

  `expr_pct` is the percentile of `expression` **within the scored batch**
  (`rank.expr_percentile`), and it buys two things a level cannot:

  - **The unit stops mattering.** TPM, FPKM and raw counts give the identical column, because a
    percentile is invariant to any monotone rescaling of abundance. A caller does not have to
    convert, and cannot convert wrongly.
  - **No imputation constant and no indicator.** A row with no expression value sits at `0.5`,
    which is what "no information" means on a percentile scale. A batch with fewer than two finite
    values — including one entirely absent — is all `0.5`: one point has no percentile.

  On the same 354,909 rows and 958 positives, at **one fewer parameter**: leave-one-screen-out mean
  0.6679 → **0.6688**, median 0.6434 → **0.6497**, twin-grouped CV 0.6434 → **0.6497**, BIC 4176.7
  → 4180.3. Against the previously shipped fit: **2 improvements, 5 ties, 0 regressions — the ship
  bar met for the first time.** Gfeller_GBM +0.0133 and IEDB_neoag +0.0165 improve; nothing
  regresses.

  **The cost, stated: the term is cohort-relative.** One peptide's `expr_pct` depends on what else
  was scored with it. `rank` and `p_response` were already per-cohort quantities; `score` is now
  one too.

  Thirteen arms were measured (`bench/results/epic_expr_arms.md`). Two negative results worth
  recording because both are the obvious idea: a **fixed constant fill does not work** — low
  TPM = 0.5 gives 0.6599 with two regressions, the GTEx typical-gene median 0.6602 with two, both
  *below* the previous fit, because a constant is not screen-neutral. And **`noncanonical`
  (missense vs other) does not replace the indicator** — a better column, present in nine screens
  of ten and 3.19% positive against 0.23% canonical, but 4,517 of 354,909 rows cannot move a
  corpus-level ranking.

### Added

- **The deposits' own gene is read before the peptide is asked.** `neoag_tested_hsa` carries
  `uniprot_id` on 99.7% of rows, `_mmu` on 97.6%, `nci_gartner_*` carry `gene_name` on 100% — and
  the pipeline read none of it, resolving the source gene by proteome search instead, which fails
  on exactly the rows that then had no expression value. The corpus builder carries both fields
  now and the accession is translated through the proteome headers. 33,197 of 363,324 rows name
  their own gene; IEDB_neoag's missing rate falls 49.8% → 43.9%, corpus-wide 5,849 → 5,774.

- **The per-allele calibration cache ships on.** A `%rank` background is a random-peptide draw
  scored under one allele's model — ~0.95 s to build, and a pure function of
  `(allele, model, background, footprint, seed, library version)`. The on-disk cache for it has
  existed, atomic and correctly fingerprinted, since it was written; it was **opt-in through
  `MHCMATCH_CALIBRATION_CACHE` and essentially nothing set it**, so every process rebuilt every
  allele it touched on every run. On the neoantigen feature build — 2,093 distinct alleles across
  fourteen workers — that was the whole cost of the stage.

  It now defaults to `$XDG_CACHE_HOME/mhcmatch/calibration`, falling back to
  `~/.cache/mhcmatch/calibration`. Set `MHCMATCH_CALIBRATION_CACHE` to share one across a SLURM
  array, or to `0`/`off`/`none`/`false` to disable. A read-only home degrades to no cache rather
  than raising. The library version is in the key, so a bump invalidates rather than serving a
  stale background, and deleting the directory is always safe.

  Measured on the neoantigen feature build, 363,324 (peptide, allele) pairs over 2,093 alleles on
  fourteen workers: the binder pass runs **1,788 s cold and 15 s warm**, a factor of 119, and the
  whole stage 2,061 s → 271 s. The two frames are **identical** — 0 of 40 columns differ over all
  363,324 rows — so this is a cache, not an approximation. Footprint after that run: 4,331 files
  and 498 MB, about 2.07 entries and 238 kB per allele.

- **`predict.binder_ranks(store, peptides, allele, …)`** — the transpose of `binder_score`: one
  allele, many peptides, which is the call shape a benchmark needs. Score-identical to
  `binder_score` by construction and pinned as such by a test. It is **not** a speed fix and is not
  described as one: measured at 1.13× on a warm allele (82,241/s against 72,966/s over 5,000
  peptides), because the cost was never the peptide loop.

- **`mhcmatch --version`.** `bench/run_epic.sh` stage 0 gates the whole reproduction on the
  installed version matching the checkout's `pyproject.toml`, and read it from this flag. There was
  no such flag, so `grep` found nothing, `pipefail` propagated, and the chain stopped at stage 0
  looking like a clean exit.

### Fixed

- **`build aggregate` has a real generator.** `bench/immuno/epic_v4_fit.py` writes the artifact into
  the library checkout directly, so the hand-copy that let the GRAND → EPIC rename reach the
  artifact but not its generator is gone. `mhcmatch build` prints the command; `bench/run_epic.sh`
  runs the chain that leads to it.
- **Two schema drifts in the artifact.** v4 moved `fit.loo` to the top level and turned `blocks`
  from a list of pairs into a dict — a consumer reading `a["fit"]["loo"]` broke, and one reading
  `[b[0] for b in a["blocks"]]` silently got the first letter of each key. Both shapes are restored,
  and `phys_scale_charge` names the second chemistry scale that is actually fitted. All three are
  now **written by `epic_v4_fit.py`** rather than patched into the file after the fact; the patched
  copy is why the shipped artifact and its generator drifted apart (below).
- **The shipped artifact is the one `bench/run_epic.sh` produces.** `epic_v4_fit.py` writes a
  *candidate*, and copying it into the library is manual, so the two can drift — and `build --check`
  cannot see it, because it compares version stamps and a hand-copied older fit stamped 0.27.0 is
  current by that test. The file this release was going to ship was fitted against a frame whose
  upstream `features_grand.parquet` and `dai_kd.parquet` were stale. Re-stamping those forced a
  rebuild that recovered an expression value for **202 of the 354,909 rows** (`expr_missing` 1.705%
  → 1.648%); every other column is bit-identical between the two frames, `mu` and `sigma` matching
  to 0 on `pres`, `occupancy`, both chemistry scales and all three corpus channels. The refit is
  better on every summary: BIC 4172.4 → 4168.6, LOSO mean 0.6602 → 0.6654, median 0.6385 → 0.6399,
  twin-grouped CV pooled AUROC 0.7665 → 0.8100.
- **The release workflow rebuilds every buildable artifact, and gates on all 27.** `publish.yml`
  ran `python tools/build_anchor_models.py` — one of the three families, through a shim — so a
  release could cut a wheel whose `corpus_tables.npz` was stamped behind `__version__` and nothing
  in the publish path would notice. It now runs `mhcmatch build` (~6.5 min: anchor ~250 s, corpus
  ~130 s) followed by `mhcmatch build --check`, which fails the release if *any* shipped artifact
  is stale, including the ones whose generator lives in the benchmark repo and which the workflow
  cannot rebuild itself.
- **`tools/build_anchor_models.py` and `tools/build_corpus_tables.py` are gone.** They were thin
  shims onto `mhcmatch._build`, and a second name for one command is a second thing to keep
  current — which is exactly how the publish workflow stayed pinned to one family. `mhcmatch build
  [target]` is the only entry point; `tools/build_recognition.py` stays, because it needs ESM2 and
  genuinely cannot run in-process. Every docstring, test message and `PROVENANCE.md` command that
  named a shim now names the CLI.

## [0.26.0] - 2026-08-22

### Fixed

- **The cassette safety screen was mis-specified for somatic neoantigens, in both of its clauses.**
  `vector.self_origin_risk` asked two questions of every unit regardless of what kind of product the
  unit encoded, and both questions were the wrong ones for a mutated product.

  *Clause 1, "target gene"* — is the unit's own gene transcribed above 0.25 TPM in an essential
  tissue — was set under **MAGE-A12 at 0.33 TPM in brain caudate**, the expression that killed two
  patients (PMID 23377668). MAGE-A12 is a **cancer-testis antigen**: a shared, *unmutated* self
  protein, so brain transcription is exactly the hazard. A somatic neoantigen encodes a sequence
  absent from normal tissue by construction and is not that object. On a 37-donor cohort the
  unconditional clause withdrew a candidate for the fact that its parent gene exists: **10 of 37
  donors lost every unit they had**, and one lost 1,098 of 1,618 to this clause alone. It is now
  gated on `Unit.kind` against the new `predict.NOVEL_PRODUCTS`; an `isoform`, a wild-type or
  overexpressed target keeps it, and **an unknown or missing kind keeps it too** — fail closed, with
  the kind carried in the reason.

  *Clause 2, "unrelated self origin"* — does a register coincide with a self peptide from an
  essential-tissue gene — read the unit's own design as the hazard. On **178 experimentally
  immunogenic somatic neoantigens** rebuilt as the 27-mer units they would enter a cassette as,
  **178 of 178 (100 %) trip it**, at a median of **36** self registers each; 36 is exactly
  `12 + 10 + 8 + 6`, the count of 8/9/10/11-mer windows of a 27-mer that cannot contain a centred
  mutation. Measured self fraction against pure geometry, per length: 60.02 / 52.6 / 44.4 / 35.2 %
  against 60.0 / 52.6 / 44.4 / 35.3 % predicted — 6,350 hits, **99.1 % of the geometric ceiling**.
  At the minimal-epitope level the clause is clean: 0 of 178 mutants in the proteome, 178 of 178
  wild types. The clause now judges only the registers that carry novel sequence, gated on the same
  `NOVEL_PRODUCTS`, with `predict.TRACT_PRODUCTS` (`frameshift`, `fusion`) novel from the variant
  offset to the end of the unit rather than at one index. `n_registers_spanning` and
  `n_hit_spanning` ride on every clause-2 reason so the exemption is auditable.

- **`vector.units_from_context` admitted only `Somatic:` windows**, discarding **317 of the 489**
  non-missense records in one real cohort and leaving the `nonconventional` quota arm nothing to
  fill itself from. All four header families are now centred: a parenthesised marker (`Somatic:`), a
  three-part `LEFT|X|RIGHT` or two-part `LEFT|RIGHT` junction (`Fusion:` / `CNV:`, truncated at a
  read-through `*`), and the trailing novel span (`Isoform:`, now parsed into a `span` field). Where
  the caller's rows carry no `kind`, the header's own product class is read rather than defaulting to
  `missense`.

- **`portfolio._ratio` raised a bare `ValueError`** naming one index. It now raises
  `portfolio.MarginalExceedsBlock` (a `ValueError` subclass) carrying the **arm**, how many of its
  units exceed `q`, and how far `--block-live` has to move.

### Added

- **`Proteome.find_exact_sources(peptides)`** — `find_sources` at `max_subs=0` without the fuzzy
  index. `_index` is a Python loop over every position of every protein (68,398,087 iterations at
  L=9, ~12.6 GB peak) and answers `<= max_subs`, which the safety screen never asks. Membership and
  provenance both come out of one sorted `(window, buffer offset)` array per length and a pair of
  `np.searchsorted` calls over the whole batch. A per-hit `bytes.find` was tried and does not scale:
  a 27-mer is native context, so 9,497 of 9,500 distinct registers hit and each would cost a full
  buffer scan.

- **`vector.screen(..., notes=[])` and a `prepare(registers)` hook on the risk callable.** The screen
  resolves the deduplicated registers of *every* unit in one query instead of one per unit
  (~19,000 calls where 1 suffices), and `self_origin_risk` caches the essential-tissue filter per
  gene instead of rebuilding it per unit and per hit gene.

- **A graded screen.** `self_origin_risk(..., graded=True, veto_tpm=5.0)` /
  `mhcmatch vector --screen-mode graded`: `min_tpm = 0.25` stays the **reporting** floor and
  `veto_tpm` becomes the **exclusion** line, so a sub-veto finding is kept as a per-unit off-target
  fingerprint rather than a refusal. `vector.offtarget_cost` turns those findings into a cost and
  `portfolio.compose(cost=..., weight_cost=...)` / `--weight-offtarget` prices it into the
  objective. **The cost is never charged to `Unit.p`**, which is a calibrated marginal `survival`
  reads literally, and `weight_cost = 0.0` is bit-identical to the previous composition. The default
  screen mode stays `veto`.

- `Composition.arms[...]` gains `mean_cost` and `max_cost`; `Composition.trace` gains `cost` and
  `cost_penalty`; the cassette report gains a `fingerprint` section for kept units.

## [0.25.0] - 2026-08-21

### Removed

- **The `pmhc = ["pandas"]` extra, which installed pandas for nothing.** It was declared as
  "convenience IO for large tables" and no module ever imported pandas for that: the only pandas in
  the package is `logo.render`, where it is **logomaker's API requirement** — logomaker takes a
  `DataFrame` — and `logomaker` already depends on it, so the `logo` extra covers it. Nothing in the
  README, docs, skill, CI or the nextflow module referenced `[pmhc]`.

  The core stays `csv` + `numpy` with **no dataframe dependency at all**, which is deliberate: a
  library should not make every `pip install mhcmatch` carry a dataframe engine to read a TSV.
  Measured on the largest table this package's own pipeline produces (22,992 rows x 43 columns,
  7.9 MB): `csv.DictReader` 0.08 s against `polars.read_csv` 0.04 s — 40 ms, inside a stage whose
  other work is 143 s. If a dataframe is ever warranted here it should be polars, but nothing
  measured so far warrants one.

### Changed

- **The shipped scorer is named `EPIC`.** Same artifact, same nine coefficients, same
  `"version": 3` — `data/aggregate_mhc1.json` now declares `"model": "EPIC"` and carries
  `"former_name": "GRAND"`, so **every recorded result under the old name is a result about this
  model** and nothing needs re-running. **E**xpression, **P**resentation, **I**mmunogenic
  **C**omplementarity names the four fitted blocks; it is not their pipeline order, which is
  presentation → expression → physchem → corpus, the two recognition blocks being the two halves of
  Complementarity.

  `rank`'s header line, the docs and the skill say `EPIC`. The model registry
  (`docs/models.rst`) keeps `GRAND` against the v2 row, which is the name that version actually
  shipped under. A consumer that asserted `components["model"] == "GRAND"` will see `EPIC`; the
  field was always display-only and no code branches on it.

## [0.24.1] - 2026-08-21

### Added

- **The corpus k-mer tables ship in the wheel.** `mimicry.corpus_counts` slid over every window of
  every reference set on every import that needed one, and for the `self` channel that is the whole
  proteome: **51.4 s** (class I, four lengths) and **14.5 s** (class II), 115.6 s for the eight
  channels the shipped model reads, **in every process**. The result is 8,000 float64s per channel
  and a pure function of the deposit — **145.4 kB of output for 115.6 s of work**.

  `src/mhcmatch/data/corpus_tables.npz` now carries all eight (class x component x species), built
  by `tools/build_corpus_tables.py` and read by `corpus_counts` for the default path in **0.002 s**.
  A custom `pmhc_dir`, `weights="locus"` or a non-default `k` each define a different table and
  still build in full — a wrong answer delivered fast is the failure mode a cache layer usually
  ships with, and there is a test for it.

  **Shipped, not cached.** A cache directory adds a staleness mode and a concurrent-write race;
  package data has neither. It is version-stamped and regenerated at release like the vendored
  anchor models, and unlike those a bump-only rebuild is genuinely bit-identical — verified against
  a live rebuild on the four deposit channels, and on totals for the four proteome ones.

### Changed

- **`load_references` builds the lengths it is asked for.** The index is per-length and so is the
  cost — ~11 s per proteome pass for `window_array`, plus **~1.0 min** at class II to resolve a
  register for each of 12,685,964 windows — and the class admits **fifteen** lengths. A run whose
  peptides are 15-mers was paying for thirteen it would never query: **~19 min against 83.1 s** for
  the one, and **65.4 s against 16.5 s** at class I. `rank --extended/--annotate` and `mimicry` now
  pass the lengths their own candidate list has. `lengths=None` still means every admitted length,
  so an existing caller is unchanged.

- **A reference source is a str list, not a padded byte array.** `_Backing` held sources as
  `dtype="S"`, which pads every row to the longest — and one thymic 9-mer (`RIHTGEKPY`, a
  zinc-finger motif) carries a `;`-joined list of ~200 accessions at **2,141 characters** against a
  mean of **7.5**. Class I alone paid **56.3 MB for 26,302 entries**, ~99 % of it NUL. Same values,
  ~2 MB.

### Fixed

- **`variant_type` was the header's provenance, not the product, so the non-conventional quota could
  never bite.** A pipeline window header carries `Somatic` in its `type` field and the consequence
  (`missense_variant`, `frameshift_variant`, ...) in `subtype`; `rank` emitted the former.
  :func:`mhcmatch.portfolio.default_arm` asks only whether a unit's kind is `"missense"`, so on any
  real donor FASTA **every** candidate was charged to `nonconventional` — the `mhc1` and `mhc2` arms
  were unfillable and `--quota mhc1=14:2,...` was unsatisfiable. Measured on a 36-donor cohort:
  5,948 of 6,437 class-I windows (92.4 %) are missense and all of them were misfiled.

  New `predict.variant_product(var)` returns the product class — `missense` / `frameshift` /
  `inframe_deletion` / ... for a `Somatic:` window, the lower-cased type (`fusion`, `isoform`,
  `cnv`) otherwise, because a fusion is non-conventional whether its junction is in frame or not.
  An unmapped consequence passes through lower-cased rather than defaulting to missense. Both
  `rank` entry points use it, and `vector.units_from_context` no longer falls back to `var["type"]`.

  Every unit test of the arms built `Unit(kind=...)` by hand, which is why all of them passed; the
  new one starts at a real header.

- **`vector --quota` composed a cassette and then built the sequence from `select` anyway.** It
  reported and did not act: `--fasta`, `--fasta-nt` and `--map` all described the n0 stopping rule's
  output, so the composition was computed and discarded. With a quota the FASTA files now carry
  **two records** — `cassette_composed` and `cassette_topk`, the same slot budgets filled by score
  alone — and the map describes the composed one. Without a quota the output is byte-identical to
  0.24.0. The report's section names are qualified (`composed:unit`, `topk:cassette`) only when
  there are two.

### Changed

- **The class-II `self` corpus table is the proteome's own k-mers, not a projected face.**
  `corpus_counts` resolved a class-II TCR face *per reference window* — `mhc2_anchors` on each one —
  and for the `self` channel that is 15 lengths x ~12.7 M proteome windows = **~192 M register
  searches**, measured at **>25 min and ~10.7 GB RSS** against **1.7 s** for the thymic deposit.
  `rank --cls mhc2 --score aggregate` was therefore not usable in practice.

  A proteome has no register, because nothing in it is presented: `thymus` and `viral` are ligand
  deposits and keep the per-window face, but for `self` the reference object is the window's own
  k-mer content. Read once at the shortest admitted length rather than through all fifteen — each
  extra length re-counts the same k-mers with a different multiplicity, which `N_k` divides straight
  back out. **14.0 s**, N = 110,932,623 windows, table fully dense.

  **Class I is bit-identical**, verified by hash on all three channels: `C_corpus_self` is a fitted
  feature and it was fitted on class-I rows, so a class-II definition may change and a class-I table
  may not.

- **`parse_variant_header` reads the three non-`Somatic` families instead of skipping them.**
  `Fusion:` and `CNV:` are colon-delimited with their own field order, `Isoform:` is pipe-delimited;
  each order was pinned against the pipeline's own `.epitopes.*.tsv` columns rather than inferred.
  They yield gene, transcript and — for `Isoform:` and `CNV:` — a real `tpm`. `Fusion:` carries
  **FFPM, not TPM**, so it lands under its own `ffpm` key and never in the slot `expr` is scored
  from. These are 7.6 % of a donor cohort's windows and *all* of its non-conventional ones, i.e.
  exactly the pool a quota's third arm draws from, and they were taking imputed expression under the
  model's largest coefficient. The returned dict is now the union of all four families' keys, so its
  shape does not depend on which header it came from.

- Nextflow: `mhcmatch_vector_quota`, `mhcmatch_vector_block_live`, `mhcmatch_vector_evenness` and
  `mhcmatch_prevalence` were read by `main.nf` and declared nowhere; they now carry defaults and
  documentation in `nextflow.config`. The `MHCMATCH_RANK` resource note still described 0.21.0's
  thymus-only corpus term.

## [0.24.0] - 2026-08-20

### Changed

- **`C_corpus` is now the exact Łuksza sum, and it costs a table lookup.** The term is
  `Z = Σ_r exp(−κ(a₀ − s(q,r)))` over a reference immunopeptidome, and mhcmatch computed it by
  walking a trie to Hamming radius 2. With an ungapped position-additive score the weight
  **factorises over positions**, so the sum over the *whole* reference set is one 20×20 matrix
  applied along each axis of a k-mer frequency table — a contraction, computed once
  (`mimicry.corpus_counts` + `mimicry.contract`), after which every query is one array index.

  | | radius-2 search | k-mer contraction |
  |---|--:|--:|
  | 340,876 queries | ~46,000 ms | **2.3 ms** |
  | `self` channel build | 75.6 s + ~7.5 GB index | **50 s**, 64 KB resident |
  | agreement with a true all-vs-all | median **0.4999** of `Z` | **5.5×10⁻¹⁶** |
  | Spearman with peptide length | **−0.502** | **+0.040** |

  The reported column is the **per-window density** `ρ = S_k/(m_k N_k) ∈ [0,1]`. `m_k` (the query's
  sliding-window count) is what removes the length artefact; `N_k` (the corpus's total window mass)
  is what puts deposits of 140,482 and 121,968,158 windows on one scale. `a₀` and the `Z/(1+Z)`
  saturation are **retired** — the first was a scale the standardizer absorbed, the second bounded
  a count that is now already bounded. `bench/results/corpus_exact.md`.

- **GRAND v3: nine terms in four hierarchical blocks.** `presentation` → `expression` → `physchem`
  → `corpus`, entered in pipeline order, so a recognition coefficient is what that term is worth
  *after* presentation and expression. Both recognition blocks are significant on entry
  (physchem LR χ²(2) = 11.0, p = 4.0×10⁻³; corpus LR χ²(3) = 15.7, p = 1.3×10⁻³) and the full model
  has the best held-out mean AUROC of any rung in the ladder (0.6927 against 0.6734 for presentation
  alone, leave-one-screen-out over nine screens). `bench/results/grand_corpus.md`.

  - `C_phys` becomes **two** columns, `C_phys_rose` and `C_phys_hydrop` (Kidera KF4). Rose carries
    the block on neoantigens; KF4 is the stronger of the two on the Chowell-family corpora that
    motivated a chemistry term in the first place. They correlate −0.836, so the statistic to read
    is the block test, not a per-term `z` — the report gives the fit in its sequential
    (Gram–Schmidt) basis and in both entry orders for exactly that reason.
  - `C_corpus` becomes **three** columns — `thymus`, `self`, `viral`. `self` and `viral` were
    dropped in 0.21.0 for what they cost (a ~7.5 GB trie), not for what they were worth; the
    contraction removes the cost. The thymus/self **sign dissociation** (+0.246 / −0.241) is intact
    and is the evidence for the biased-sample mechanism.
  - `C_corpus_missing` is **removed**. It flagged a peptide with no cache entry; the exact sum has a
    value for every canonical peptide of an admitted length, so the column would be identically
    zero.

- **`complement.burial` averages over the TCR face instead of summing.** The face is `L − 5`
  residues wide and the Rose scale is strictly positive (0.52–0.91), so the summed column was
  **Pearson +0.954 with peptide length** (n = 60,000) — a chemistry term that was 91 % ruler. The
  averaged column sits at −0.010, its marginal within-screen AUROC rises 0.5098 → 0.5646, and the
  two scales become comparable for the first time (their correlation moves from −0.20 to −0.836,
  which is what they always were under the length variance). `burial(..., per_residue=False)`
  reproduces a pre-0.24.0 number.

- **`mimicry.SHAPES` is one `κ` per component, not a `(κ, a₀)` pair**, profiled by
  `bench/immuno/corpus_exact.py`: `thymus` 3.0, `self` 5.0, `viral` 8.0. `mimicry.RADIUS` and
  `mimicry.corpus_radius()` are removed with the search. Writing `γ = (1−e^−κ)/(1+19e^−κ)`, the
  contraction keeps an order-`j` interaction at exactly `γ^j` — so `κ` is a single scalar bandwidth
  running from pure composition to exact k-mer matching, and the three components sit at
  γ = 0.49 (graded tolerance), 0.88 (unidentified) and 0.99 (near-exact). Derived and verified to
  10⁻¹⁵ in `bench/results/kmer_spectrum.md`.

- **`rank` is the rank by score**, dense and 1-based, rather than the row's position in the file.
  Known epitopes are still floated to the top of the listing; that is a display choice and no longer
  renumbers the ranking.

### Added

- **`p_response`: the score on a probability axis.** `rank.probability` picks the single additive
  offset that makes the mean fitted probability over the pool equal a declared prevalence, and
  reports `σ(s + b)`. `--prevalence` (default 0.0602 = TESLA's 37 of 615). It is a **prior shift,
  not a recalibration**: additive in log-odds, so it moves no rank. What it buys is portability —
  a raw-score cut-off means nothing across cohorts whose base rates differ by four orders of
  magnitude, and "P ≥ 0.2 at an assumed 6 % pool prevalence" is a statement another cohort can be
  held to.

- **`portfolio.compose`: cassette composition to quotas.** `{"mhc1": (8, 2), "mhc2": (4, 1),
  "nonconventional": (3, 1)}` reads *eight class-I slots, of which at least two should respond*.
  `mhcmatch vector --quota`, which reports the composed cassette **and the same slot budgets filled
  by score alone**, side by side. `P(≥ k)` is not a modular set function whenever two units share a
  block, so no pointwise score can be sorted to maximise it: on nine candidates at target 1, the
  composer reaches P(≥ 1) **0.6550 against top-4-by-score's 0.4806**, and spreads over four
  allotypes instead of one without being told to. At target ≥ 2 it *concentrates*, and is right to —
  two units in two blocks need both blocks live.

  Arms are disjoint by construction (`Unit.kind`): a frameshift product is charged to
  `nonconventional` rather than to `mhc1`, or "at least one non-conventional epitope responds" would
  be satisfied for free and the quota would never change a cassette.

- **`portfolio.survival` and `portfolio.coverage`.** `survival` is the whole tail of the block
  response model, **exact** by convolution over blocks — `O(B m²)`, no `2^B` live-set enumeration
  and no Monte Carlo. `p_at_least` reads it and has lost its `n_mc`/`seed` arguments. `coverage`
  gives the Gini index and H/H_max of allotype spread, over the **donor's own distinct allotypes**:
  a patient homozygous at *B* has five, and scoring an even cassette against a denominator of six
  would report a genotype as a design flaw.

- **`variant_type`** on the `rank` output, from the FASTA header's `type` field or the table's
  `type` column. Reported, never scored — it is what lets the cassette layer hold a quota of
  non-conventional epitopes.

- `mimicry.corpus_counts`, `mimicry.contract`, `mimicry.face_kmers`, `complement.PHYS_SCALE_HYDROP`,
  `rank.AGGREGATE_BLOCKS`, `rank.PHYS_COLUMNS`, `rank.POOL_PREVALENCE`.

### Removed

- **`$MHCMATCH_REFERENCE_CACHE`, and `load_references(cache=)`.** There is nothing left to cache on
  the corpus path: the thymic table builds in 0.5 s and the contraction in ~1 ms. **If a cluster has
  been pointing this at shared storage, that directory (~1 GB) is now dead and can be deleted** —
  the variable is simply unread. `mimicry.CACHE_VERSION` and the whole cache read/write path go with
  it (~130 lines). The indexed search remains for `features()`, `annotate()` and the self-mimicry
  safety scan, which report *which* reference was hit and cannot be answered by a weighted sum.
  Removed from `slurm.config` and the nextflow README.

- `mimicry.RADIUS`, `mimicry.corpus_radius`, and the artifact's `corpus_radius` key.

### Note

- The bump invalidates `$MHCMATCH_CALIBRATION_CACHE` entries (`predict._fingerprint` carries
  `__version__`). That cache is safe to leave in place; it rebuilds.

## [0.23.0] - 2026-08-20

### Added

- **`--core`: the binding core, on every output that can carry one.** `rank`, `predict` and `neoag`
  gain the flag; the cassette map (`vector --map`) carries `core` unconditionally, beside the
  `core_start` / `core_end` it already had but never had residues for. Three columns — `core`,
  `core_offset`, `core_source` — following NetMHCpan's definition, "the minimal 9 amino acid binding
  core directly in contact with the MHC (i.e. excluding potential insertions)", with `core_offset`
  its `Of`, 0-based.

  **The core is residues, never a padded frame.** The parenthesis in that definition is the operative
  part: where an alignment to a 9-mer motif needs a gap, the inserted position is not part of the
  core. So it is nine residues whenever the peptide can fill nine — every class-II core, and a
  class-I 9-, 10- or 11-mer — and the peptide's own residues when it cannot. A gap character would
  not be neutral in an amino-acid column in any case: `B` is Asx in IUPAC, so a reader would take it
  for a real ambiguity code.

  **Class I holds both anchors and lets the middle give way.** `mhcmatch.store.binding_core` resolves
  `diffusion.MHC1_CORE` through `store.mhc1_positions` — the same mapping the scorer uses, so the
  reported core is the residues the model actually read. A 9-mer is its own core; a 10- or 11-mer
  drops one or two central residues (NetMHCpan's `Gp`/`Gl` deletion); below nine the `+5` and `-4`
  positions collide and the losing *slot* is dropped rather than padded — not a residue — so every
  residue still appears exactly once and an 8-mer's core is the 8-mer.

  **Class II is the register-anchored 9-mer**, matching NetMHCIIpan's `Core`/`Of`. Which register
  produced it is a column rather than a footnote, because the two disagree often on real ligands:
  `core_source` reads `model` where the per-allele `AnchorModel.best_register` was used (`predict`
  and `rank fasta`, which already computed it to score with and until now threw it away),
  `heuristic` where there is no allele (`neoag`, `rank table`), and `footprint` for class I, where
  there is no register to choose. A core nobody can attribute is not auditable.

  Free — nothing new is computed — reported, and **never scored**: the aggregate reads the peptide,
  and `--core` cannot move a ranking. `predict --core` is distinct from `--footprint core`, which
  changes what the model scores; the help text says so.

  Nextflow: `params.mhcmatch_{predict,rank,neoag}_core`. The `MHCMATCH_RANK` stub now passes `core=`
  to `rank.columns()` rather than hardcoding a header, which is the module's own rule.

### Changed

- `write_scored_csv` takes `core=`. **The 57-column `.epitopes.scored.csv` schema is unchanged by
  default** — it is a contract with the downstream pipeline modules, and `DictWriter`'s
  `extrasaction="ignore"` would have dropped a stray key silently rather than failing, so widening
  it is the caller's explicit choice and nothing else.

## [0.22.0] - 2026-08-20

### Removed

- **`mhcmatch.ipred` is gone, and so is `mhcmatch/data/ipred_mhc1.json`.** The legacy
  physicochemical immunogenicity predictor — three features summed over the whole peptide (`pc1`,
  `pc2`, `length`), two class-conditional diagonal Gaussians fitted by EM, a two-parameter Platt
  calibration, **13 fitted parameters** — shipped from **v0.9.0** and is removed here. **0.21.0 is
  the last released version that carries it.** `from mhcmatch import ipred` now raises
  `ImportError`; the whole `__all__` (`PARAMS`, `feature_names`, `features`, `score`, `log_p`,
  `p_immunogenic`, `residue_scores`, `parameters`) and the `python -m mhcmatch.ipred` entry point
  go with it.

  It was never a term of any shipped scorer, and it does not earn one. On **355,052 peptide ×
  allele rows / 1,101 immunogenic over 10 screens** of the cleaned grand corpus, one unpenalised
  intercept per screen (`bench/results/neoag_ipred_vs_complement.md`):

  | model | `ipred` z | within-screen median AUROC | BIC |
  |---|--:|--:|--:|
  | `BOECRT` (shipped) | — | **0.6504** | **4201.7** |
  | `BOECRT` + `ipred` | +0.22 | 0.6506 | 4214.5 |
  | `ipred` instead of `complement` | +1.12 | 0.6399 | 4218.4 |

  Adding it moves within-screen median AUROC **+0.0002** (0.6504 → 0.6506) and worsens BIC
  **+12.7** (4201.7 → 4214.5), against a per-column penalty of log(355,052) = 12.78. Swapping it in
  for complementarity costs **0.0105** AUROC (0.6504 → 0.6399) at z = +1.12. It is not redundant —
  r(`ipred` log-odds, `complement` log-odds) = **+0.2018** over those 355,052 rows, and **+0.2045**
  over 362,324 human rows in a separate measurement, with all six `complement` blocks together
  explaining R² = **0.5113** of its variance (`bench/results/ipred_residual.md`) — the variance it
  carries alone simply does not help once complementarity is present.

  **The record survives in full**, in `docs/complementarity.rst`, section "`ipred`: the retired
  predecessor": the shipped configuration's pooled out-of-fold AUROC **0.712** / macro **0.607** on
  peptide-grouped 5-fold CV over **694,507 rows / 35,595 immunogenic across 7 label sources**; its
  **+0.059** pooled / **+0.030** macro margin over the summed-Kidera baseline (Chowell 2015 /
  Pogorelyy 2018) at 0.653 / 0.577, and that baseline's own win on `chowell` at **0.703** vs 0.680
  on 9,806 peptides / 5,035 immunogenic; calibration Brier **0.2282**, ECE **0.0661**, AUROC
  **0.6804** on n = 9,806; parameter stability `d[pc1]` **+0.2395** CI [+0.1844, +0.3020] and
  `d[length]` **−0.2189** CI [−0.2692, −0.1632] over 1,000 peptide bootstraps; human↔mouse transfer
  **0.733**/**0.648** human-trained vs 0.654/0.625 mouse-trained on 649,466 / 45,041 rows; and the
  two cohorts where it **beat** complementarity as a single unfitted feature — VACCIMEL AUROC
  **0.6324** vs 0.5774 (93 rows / 27 immunogenic) and GBM **0.6450** vs 0.6186 (109 / 26).

- **`rank` output loses the `physchem_ipred` column.** `rank.BASE_COLUMNS` drops it, so
  `rank.columns()` returns one fewer column in every mode and every `mhcmatch rank` TSV is one
  column narrower; the `Ranked.physchem_ipred` field is gone, so `Ranked(..., physchem_ipred=...)`
  is now a `TypeError`. `Ranked.physchem` — the complementarity recognition term — is a different
  field and stays. The Nextflow `MHCMATCH_RANK` stub derives its header from `rank.columns()` at
  runtime and needed no edit, which is what that design was for.

- **`mhcmatch explain` loses two outputs.** `explain --peptides` drops the `ipred_logp` TSV column
  (13 columns → 12); single-peptide `explain` drops the printed `ipred log P` line. No flag or
  subcommand is removed.

### Changed

- **The property basis outlived the artifact that first carried it.** `PROPERTY_PC1` /
  `PROPERTY_PC2` in `mhcmatch/data/aa_tables.py` are now sourced, in code and in
  `data/PROVENANCE.md`, to what they actually are: **derived/computed, label-free** — the first two
  principal components of the 20 × 142 property matrix, column-standardized over the 20 residues,
  by SVD, regenerated with `python bench/ipred/pca.py`. PC1 carries **32.79 %** of total variance
  with residue order `I F L W V M C Y A P G T H S Q N E K D R`, PC1+PC2 **51.2 %**, 10 components
  **91.3 %** (`bench/results/ipred_pca.md`). They previously pointed at `ipred_mhc1.json`
  (`residue_scores`), a file that no longer ships.

- **`data/PROVENANCE.md` keeps the `ipred_mhc1.json` entry, marked retired.** A result recorded
  against 0.21.0 or earlier cites a file this package used to carry; a provenance record that
  deletes retired artifacts cannot say where that number came from.

- **The letter `V` stays defined.** `BDEVF` is a published model name with recorded coefficients —
  `ipred` at **+0.2707**, 95 % bootstrap CI [+0.2349, +0.3075], second of seven, on 16,802
  peptide × allele rows / 8,258 immunogenic (`bench/results/neoag_glm.md`) — and `mhcmatch.mimicry`
  is documented as fitted residual to it. `V` was always named after the *generation* (vanilla
  physicochemistry), not the module, which is exactly what lets the name survive the removal.

- **The shipped model is untouched.** `GRAND`'s seven terms, its coefficients and
  `data/aggregate_mhc1.json` are unchanged; `physchem_ipred` was a reported column, never a
  feature. Vendored anchor models are re-stamped for the version bump.

## [0.21.0] - 2026-08-20

### Changed

- **The shipped model is `GRAND`: seven terms, and Complementarity is exactly two factors.**
  `rank` scored `BOECRT` while the benchmark had moved on. Fitted by
  `bench/immuno/grand_corpus.py` over 354,909 rows / 958 positives across nine screens, one
  unpenalised intercept per screen and no global one, BIC 4160.1:

  | term | coefficient | z | ΔBIC if dropped |
  |---|--:|--:|--:|
  | `expr` | **+0.3250** | +5.99 | +23.3 |
  | `C_phys` | +0.2579 | +4.30 | +5.7 |
  | `C_corpus_thymus` | +0.1871 | +5.48 | +6.1 |
  | `binder` | +0.1408 | +4.01 | +3.5 |
  | `occupancy` | +0.1072 | +5.33 | +14.6 |
  | `expr_missing` | +0.0994 | **+6.37** | **+28.1** |
  | `C_corpus_missing` | −0.3510 | −3.96 | +7.0 |

  Leave-one-screen-out median AUROC **0.6391** (0.5174 VACCIMEL to 0.8744 Neopep), scored with the
  mean intercept — what a new cohort gets in deployment. Coefficient signs stable over 400
  **patient-cluster** resamples (rows from one patient share a tumour, an HLA type and a
  sequencing run, so a row bootstrap reports intervals several times too narrow).

  The four recognition columns collapse to two, and **neither is fitted on neoantigen labels**:

  - **`C_phys`** — `complement.burial`, the Rose burial propensity summed over the TCR face.
    Chosen from 576 candidate scales by ΔBIC *inside this model*, not by standalone AUROC. An
    imported basis means **zero fitted residue parameters**, and it carries a cysteine loading of
    **+0.108** against the retired `complement`'s **+0.693**.
  - **`C_corpus_thymus`** — `mimicry.corpus_R` on the thymic channel. Positive, because the thymic
    immunopeptidome is a *biased* sample of self: mTECs promiscuously express tissue-restricted
    antigens under *Aire* and *Fezf2* precisely to purge the clones worth purging. Thymic ligands
    score 0.7303 against 0.7222 for non-thymic **presented** self (Cohen's *d* = +0.1650,
    p = 1.0×10⁻⁸⁰, presentation held constant). The `thymus`(+)/`self`(−) sign dissociation is the
    evidence for the mechanism; no single-mechanism account produces it.

  Every alternative was measured and each costs BIC to add back: Kidera KF4 +9.0, KF2 +12.8, the
  `self` corpus channel +8.1, `viral` +11.6, `viral_R` +11.6, `C_aa` +6.7.

- **An aggregate score no longer needs the host-proteome index.** `self_tcr` was `BOECRT`'s
  second-largest coefficient, so scoring forced 6 min 15 s and ~7.5 GB — the largest single cost in
  the package — and `--no-self --score aggregate` had to be refused outright. `GRAND`'s corpus term
  is thymic only (26,513 peptides). Measured uncached: **2.0 s** to build the index, **0.22 ms per
  peptide** to query, **11.8 s** for a whole `rank table` run. The refusal is gone and the Nextflow
  `MHCMATCH_RANK` profile now sizes on `--extended`/`--annotate`, which is what actually loads the
  self reference.

- **An unmeasured component reports NaN, not 0.** Under `allow_missing` a component with no
  reference index standardized to the training mean and contributed exactly `coef × 0` — printing
  as `0` on `autoimmune`, which reads as "no self-similarity found" when the truth is "never
  looked". `logodds` and `autoimmune` are NaN when any component is absent, because the fitted
  coefficients describe the full set.

### Fixed

- **`C_corpus` is the Łuksza form with its length term.** `mimicry.corpus_R` shipped
  `Z = Σ_d n_d e^{+k d}` against the column the model was fitted on,
  `Z = Σ_d n_d e^{-k(a_0-(L-d))}`. Two defects in one expression: the sign on `d` was flipped, so a
  neighbour 2 substitutions away weighed e^{2k} = 90× *more* than an identical one instead of 90×
  less; and the `e^{k(L-a_0)}` factor was dropped on the reasoning that `a_0` is unidentified. It
  is — *at fixed length*. Peptide length varies across a real corpus, so that factor is a genuine
  per-row term spanning the same 90× between a 9-mer and an 11-mer.

  Measured over the 328,276 peptides with a thymus `tcr5` cache entry, the shipped variant was a
  different column and not a rescaling: mean *R* 0.771 with 77.2 % of peptides above 0.5, against
  the fitted mean of 3.29×10⁻⁵ with none above 0.5 (fitted *Z* stays under 1.320×10⁻³, the linear
  regime). Spearman +0.705, Pearson +0.448. Fitted shapes now vendored as `mimicry.SHAPES`.

- **`mimicry.masks` is class-aware, so `cls` stops being a parameter that does nothing.** `masks`
  took only a length and hardcoded the class-I anchor set; `corpus_R` accepted `cls` and never
  passed it anywhere. Every class-II ligand was read on the class-I layout — for
  `AAAKFVAAWTLKAAA` the real anchor set is `[4, 7, 9, 12]` against the class-I `[0, 1, 2, 13, 14]`,
  sharing **no position at all**. `masks(length, cls, peptide, register)` delegates class II to
  `complement.mhc2_anchors` and refuses to guess a register. `complement.burial` already passed
  `cls` through, so `C_phys` was correct at class II throughout.

- Bytecode under `src/mhcmatch/data/__pycache__/` is no longer tracked. The `!src/mhcmatch/data/**`
  negation — which exists so a model file named `.fasta` cannot be swallowed by the `*.fasta` glob
  — was re-adding five `.pyc` files to the index, producing a spurious diff on every interpreter
  version. Both halves are now pinned by a test. hatchling already kept them out of the wheel, so
  this was index hygiene, not a release defect.

### Added

- `complement.burial(peptides, cls=, scale=, registers=)` — `C_phys`. `scale=` reaches the Kidera
  factors (`"KIDERA:KF4"`, `"KIDERA:KF2"`) for comparison against the Chowell-family literature;
  they lose to `"Rose"` on the neoantigen corpus, and a number produced with a different scale is a
  **comparison**, never a silent substitution.
- `mimicry.corpus_R(peptides, refs, cls=, shapes=, radius=, components=, registers=)` — `C_corpus`.
  `components=` reaches the `self` and `viral` channels for the ladder; only `thymus` earns its
  parameters inside the model. Report the ladder anyway — the sign dissociation is the evidence.
- `mimicry.SHAPES`, `mimicry.RADIUS`, `mimicry.corpus_shapes()`, `mimicry.corpus_radius()`,
  `mimicry.NEOAG_COLUMNS`, `luksza.SHAPE`, `rank.CHANNEL_COLUMNS`.
- **`mhcmatch.portfolio` — cassette composition read as a portfolio rather than a ranking.** It
  computes nothing new about a peptide: it takes the scores the rest of the library produces and
  says what a proposed *set* of them is worth. Ten functions, none of them fitted —
  `pareto_front`, `nondominated_rank`, `crowding_distance`, `linearly_supported`,
  `chebyshev_score`, `corner`, `p_at_least`, `n_effective`, `dispersion`, `betabinom_rho`.
  `vector.select` stays the rule; this is the diagnostics. `docs/portfolio.rst` and
  `notebooks/09_cassette_composition.py`.
- `docs/burial.rst` and `docs/corpus.rst`, with the `C_corpus` formula, its five steps and the
  fitted shape table.

## [0.20.0] - 2026-08-20

### Changed

- **A model now reports the features it used, and refuses to run without them.** `BOECRT` declares
  nine features. Four of them — `viral_R`, `viral_tcr`, `self_tcr`, `thymus_tcr` — were **never
  written by anything in the package**, so `aggregate_score` substituted their training means and
  each contributed `coef × 0` to every candidate. `--extended` did not repair it: the CLI computed
  those channels at `cli.py:465`, *after* `rank_fasta` had already scored, and only printed them.
  **`mhcmatch rank` therefore scored `BOEC` on every run, with or without the flag, while reporting
  `BOECRT`** — and the Nextflow `RANK`/`VECTOR` path inherited it.

  Three measures of the damage, which disagree, and all three are worth knowing:
  **38.0 %** of the model's total absolute weight was inert (`sum |coef| = 1.3875`, of which
  `self_tcr` alone is +0.3154 — its second-largest coefficient, above `complement` at +0.1790 and
  `binder` at +0.1418); the emitted **ordering was unaffected**, because a constant offset cannot
  reorder anything, so no shipped output was numerically wrong; and the **accuracy cost is +0.008
  AUROC** (`BEC` 0.6628 → `BECRT` 0.6707, within-screen median over 7 screens,
  `bench/results/neoag_cohorts.md`). What was wrong was the reported model, not the ranking.

  `rank_fasta` / `rank_table` take a `channels` callable that supplies the four before scoring;
  `aggregate_score` raises, naming the feature, when a declared one is absent.

- **`rank` no longer falls back to the gate in silence.** The whole aggregate branch sat inside a
  bare `except Exception: score = "gate"`, so a missing artifact or an absent numpy swapped in a
  different model — a two-term noisy-AND returning a probability where the aggregate returns
  log-odds — said nothing, and left `components["model"]` unset. Asking for the aggregate and
  getting the gate is not a degraded answer; it is a different one.

- **`--no-self` cannot be combined with `--score aggregate`.** The self-mimicry reference supplies
  `self_tcr`. The combination is refused before any work starts, naming the feature. Use
  `--score gate`, which does not use it.

- **`rank --score aggregate` now costs what the model costs** — but far less than it used to, and
  once rather than every run. Building the self-proteome reference was 6 min 15 s; it is now
  **75.6 s**, and **0.82 s** from cache (92×). Before 0.20.0 it was free only because four of the
  nine features were never computed. The Nextflow `MHCMATCH_RANK` process is sized accordingly
  (16 GB / 4 h, or 8 GB / 1 h under `--score gate`).

- **A non-finite value inside a supplied column is imputed *visibly*, not silently.** One candidate
  with no IC50, or a frameshift with no wild type, still takes the training mean — dropping it would
  lose a real candidate — but the new `imputed` column names which features that row had to impute.
  A placeholder nobody can see is the defect; a placeholder the row declares is incomplete data.

- **`rank.columns()` takes `score=`.** The header carries the model's own features: the four
  recognition channels are columns when the aggregate scored and absent when the gate did.
  `binder` — a model feature — was missing from the header entirely and is now in it.

### Added

- **A reference cache, `$MHCMATCH_REFERENCE_CACHE`.** Point it at a directory and the built mimicry
  indexes are written once and memory-mapped thereafter: **0.82 s against a 75.6 s build, 92×**.
  Point it at *shared* storage and a Nextflow or SLURM fleet builds once and every task loads in
  under a second — and tasks co-resident on a node share the mapped pages through the OS page cache
  instead of each holding its own ~7.5 GB copy. Entries are keyed on the reference files, the
  channel projection and `CACHE_VERSION`, so a changed input rebuilds rather than being trusted.
  `seqtree.Index` already had `save`/`load`; the representative windows and sources are two
  `|S` arrays loaded with `mmap_mode="r"`, and `features()` touches them only for a best hit, so
  the access stays sparse.

- **`Proteome.window_array(L)`** — the vectorized counterpart of `windows(L)`: a sorted `|S{L}`
  array via `sliding_window_view` + one `np.unique`, replacing a per-window Python loop that ran
  `all(c in _AA for c in w)` over 12 M windows × 9 residues. **11.0 s against 30.0 s** on the human
  proteome for 12,073,995 distinct 9-mers, identical output. `windows()` still returns a `set` for
  the O(1) membership its own callers need.

  Packing residues into `uint64` (5 bits each) and sorting integers instead was tried and is **4×
  slower** — 44.5 s — because the shift/or loop costs more than numpy's fixed-width byte sort saves.
  Measured, not assumed.

- The per-channel projection in `load_references` is now one `np.unique` over a fixed-width byte
  view rather than a `setdefault` over 12 M strings. Because the window array is sorted, the first
  occurrence of a projection is the lexicographically smallest window carrying it — the same
  representative the old loop chose, so the features are byte-identical.

- **`n_alleles_presenting` / `alleles_presenting`**: how many of the queried allotypes present this
  peptide, and which, banded on the presentation %rank at `--rank-threshold`. `predict_windows`
  already scored every allele and kept only the best, so this is free. A peptide presented by three
  of a donor's six class-I allotypes is a different bet from one presented by one — in the response
  model of `mhcmatch.portfolio` it spans three blocks by itself. **Column only**, not a fitted term.

- **`physchem_ipred`**: `mhcmatch.ipred.log_p` as a reported column, explicitly **not in the
  model**. It is the best single feature on both cohorts where the fitted aggregate sits at chance
  (VACCIMEL AUROC 0.6324 on 93 rows / 27 positives; GBM 0.6450 on 109 / 26 — against `binder` at
  0.5065 and 0.5767; `bench/results/neoag_cohort_scan.md`), which is worth being able to see.

### Fixed

- `aggregate_score`'s docstring quoted `viral_R`'s fitted sigma as `3.8e-8`; the shipped artifact
  has `4.729e-11`, three orders of magnitude apart. The docstring was stale.

## [0.19.0] - 2026-08-19

### Changed

- **`mhcmatch rank` now scores with the fitted model.** `rank._finish()` had always scored with the
  two-term noisy-AND `gate_probability(presentation, physchem)` while the fitted aggregate sat in
  `data/aggregate_mhc1.json` with **zero internal callers** — so the shipped ranking and the
  published coefficients were two different models. The aggregate is now the default;
  `--score gate` reproduces the old ordering for comparison.

- **The vendored aggregate is `BOECRT`,** refitted on the cleaned grand corpus (see
  `bench/results/neoag_aggregate_boecrt.md`). `O` (occupancy) replaces `D` (agretopicity): the
  ratio does not resolve in any of seven parameterisations, and the one that appears to is 0.9955
  correlated with mutant affinity. **Occupancy is the strongest term in the model at z +5.16**,
  above `binder`'s +4.00 — the two share the affinity axis and the fit splits it between them, with
  occupancy taking the larger share.

  355,052 rows / 1,101 positive / 10 screens, within-screen median AUROC 0.6504. The corpus is the
  cleaned one — pathogen epitopes and unmutated self windows removed, host keyed on the MHC genus,
  label conflicts counted, CEDAR and Gfeller held out — rather than the uncleaned set the previous
  artifact used, and the Łuksza shape is the refitted `tcr5`/k=2.25/a₀=20.0 rather than the
  hardcoded 1.0/24.0. **`BOECRT` is therefore not comparable to the older `BECRT` record (0.6707 on
  7 screens): different corpus, different population.** The generator refuses to write if a Sahin
  peptide reaches the fit; it checked 8 and found 0.

- **`Ranked.binder` is written explicitly.** `rank_fasta` set `presentation` from the presentation
  head while `rank_table` set the same field from the *binder* rank — two different quantities under
  one name, which the aggregate would have read as the wrong feature. Both now write `presentation`
  and `binder` separately.

### Added

- **`docs/neoantigen.rst`** — the scorer documented end to end: every term, what it was fitted on,
  why occupancy replaces agretopicity, how to read the output, and a limits section that states the
  held-out numbers and the two mimicry terms whose direction is not established.
- `mhcmatch rank --score {aggregate,gate}`.

## [0.18.0] - 2026-08-19

### Added

- **`rank.occupancy` and the `occupancy` column — agretopicity, taken from the binding equilibrium
  instead of from a ratio.** `mhcmatch rank` now emits the fraction of MHC a peptide holds,
  `a/(1+a)` with `a = [P]/Kd` and `[P]` = `rank.PEPTIDE_NM` (10 nM).

  The benchmark's fitted aggregate carried `dai` = `log10(Kd_WT/Kd_MT)` with a **negative**
  coefficient and an interval crossing zero, and the marginal was negative too — within-screen
  median AUROC 0.4986 against the binder %rank's 0.6383, in nearly every screen. The cause is that
  the raw ratio is not what Łuksza et al. (Nature 2017;551:517, doi:10.1038/nature24473) use: they
  apply a pseudocount ε = 1/3687 nM to both dissociation constants and exclude substitutions away
  from P2/PΩ. `AffinityModel.amplitude()` has shipped that pseudocount since 0.9.0 (eq. 9) and no
  fit had ever called it.

  Restoring both corrections is not enough — seven parameterisations were fitted and none resolves.
  The ε sweep appears to (z climbs to +3.32) but the term is then 0.9955 correlated with
  `-log10(Kd_MT)`: it improves by *deleting* the agretopicity, leaving mutant affinity under another
  name. The competitive Langmuir occupancy is the term the equilibrium actually supplies, and three
  things we had been bolting on by hand fall out of it — the binder gate (a non-binding mutant
  occupies nothing whatever its wild type does), the pseudocount (the free-MHC `1` **is** Łuksza's
  ε, which is why that ε has units of inverse concentration), and a bounded scale whose steepness is
  fixed at 1 rather than fitted.

  Occupancy is **absolute** where the binder %rank is allele-relative, so the two are additive, not
  redundant: fitted together, `binder` holds z +6.5 and occupancy carries z +3.6 to +3.8, stable
  across `[P]` from 1 to 1,000 nM — it is not fitting its own concentration parameter. Model
  within-screen median AUROC 0.6526 against 0.6390. And it needs no wild type, so it is defined for
  a frameshift or fusion product, where `dai` had to be fabricated from a per-cohort q90 quantile.

### Changed

- **`agretopicity` is reported, not fitted.** It stays in the `rank` output and is no longer a term
  of the aggregate. `MODELS.md` names the model `BOECRT`.
- **`rank.BASE_COLUMNS` is 13 columns**, `occupancy` inserted before `agretopicity`. The Nextflow
  stub reads the header from `rank.columns()` so it tracks the change without editing.

## [0.17.0] - 2026-08-19

### Added

- **`mhcmatch.luksza` — the `R = Z/(1+Z)` recognition term, so `viral_R` is computable in-library.**
  The fitted aggregate carries a `viral_R` coefficient, but the Boltzmann sum lived only in the
  benchmark repo: `rank.aggregate_score` was a public function with a feature no installed user
  could supply. `luksza.viral_r(peptides)` now goes from peptides to that column end to end against
  the same viral ligandome and radius the coefficient was fitted with.

  `k` and `a0` are **read from the shipped artifact**, not hardcoded, so a refit needs no code
  change. `r_term` reproduces the benchmark's implementation **bit-identically** (0.0 over 200
  random trials, asserted in the suite), which matters because the coefficient is only meaningful
  against the exact sum it was fitted on.

  Sanity, on the shipped reference: the two genuine viral epitopes among five test peptides —
  GILGFVFTL (influenza M1) and SLYNTVATL (HIV gag) — score 3.9e-07 and 1.5e-06 against 2.1e-08 to
  3.2e-08 for the tumour peptides, an order of magnitude apart in the right direction.

  **Speed was measured, not assumed.** End to end on 20,000 peptides against the 57,331-peptide
  viral set at radius 4: 57,000 peptides/s, of which the seqtree neighbour search is 98.6 %,
  `counts_by_distance` 1.3 % and `r_term` 0.1 %. Both new functions are already an order of
  magnitude faster than the search feeding them, so they are deliberately left un-vectorised and the
  module says so — optimising them would buy 1.4 % of the run.

- **`docs/cli.rst`** — the CLI had no reference page anywhere in the docs, only a bare
  comma-separated list in the skill and prose scattered through the README. Nineteen commands
  grouped by task and by axis, the batch/`--threads` rule, the `predict` vs `restriction`
  distinction, and the environment variables a cluster needs. Wired into the toctree with its own
  landing card, alongside a card for the safety page which also had none.

- **The Nextflow `MHCMATCH_VECTOR` process emits the cassette map** (`map` / `map_json` channels,
  on by default, `optional: true`). 0.16.0 shipped the map and the module could not produce it, so
  the pipeline was a release behind its own library.

### Measured

- **Regenerating the vendored anchor models for a version bump moves no prediction.** The load guard
  keys on `__version__`, so every release rebuilds them; 0.16.0 -> 0.17.0 leaves `panel_sha` and
  `params` unchanged and the refit is deterministic, giving **bit-identical scores** (max
  |old - new| = 0.0 over 12,000 scorings across four alleles). The guard is a provenance guard, not
  a correctness one, and a downstream deliverable does not need re-running for a release that
  changes neither the panel nor the parameters.

### Fixed

- **`recognition.default_head` returns `complement`, not `posbayes`.** The README and
  `complementarity.rst` both stated the old value; `lowest_bic_head` is the one that still reports
  `posbayes`, and the two answer different questions, so both are now shown together. Caught by a
  clean-install smoke test against the published 0.16.0 wheel.

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

- **Notebook 8 — `08_ranking_and_cassette.py`**, the applied pipeline end to end on a **mock scored
  table**: `rank_table` recomputing presentation rather than re-sorting an upstream column, why a
  cassette unit is the 27-mer window and not `rank`'s 9-mer, `screen` withdrawing before capacity is
  spent, `select`'s explicit `n0`, `order` choosing no spacer at all, and `epitope_map` marking the
  three junction-spanning epitopes that belong to no gene. Runs in ~18 s with no cohort data.

- **`docs/cli.rst`** — the CLI had no reference page anywhere in the docs, only a bare
  comma-separated list in the skill and prose scattered through the README. Nineteen commands
  grouped by task and by axis, the batch/`--threads` rule, the `predict` vs `restriction`
  distinction, and the environment variables a cluster needs. Wired into the toctree with its own
  landing card, alongside a card for the safety page which also had none.

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
