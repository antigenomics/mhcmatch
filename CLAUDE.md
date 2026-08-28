# CLAUDE.md — working conventions for mhcmatch

**Authoritative context lives elsewhere — read it first:**
- [`ROADMAP.md`](ROADMAP.md) — the agent contract: what mhcmatch is, phase status, open loops.
- `../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` — the method/statistics spec (manuscript repo).
- `../../manuscripts/2026-mhcmatch/results/CURRENT.md` — every published number and which artifact
  it belongs to; `issues_major.md` beside it holds the decisions that are the author's, and
  `issues_minor.md` the repairs already applied.

**The shipped scorer is artifact version 10 in library 1.4.0** — nine fitted terms, `binder` as the
fitted presentation term (`issues_major.md` M1 and M12, both taken by the author on 2026-08-25). Its
own `verdict` block reads `"ship": false`, which is worth knowing before quoting it: v10 was shipped
on the author's word on 2026-08-28 *against* that bar, over two thin regressions (IEDB_neoag and
ITSNdb, the latter near chance under v9 too). Do not replace
`src/mhcmatch/data/aggregate_mhc1.json` without the author's word: it moves every number in the
manuscript, and `build --check` cannot see that it changed — but
`test_the_shipped_artifact_is_pinned_to_the_fit_that_produced_it` now can, by digesting
`(coef, mu, sigma)`.

**A pooled null must leave the queried allele out.** `background="ligand"` pooled the residue
marginal over *every* allele's ligands, the queried one included -- a rounding error on a balanced
panel and the entire score on a skewed one. `H-2-IAb` is **6,483 of 6,705** mouse class-II ligands
(96.7%), so its null was its own motif and the benchmark reported frequent AUROC **0.322**, below
chance on 6,483 ligands. Fixed in 1.5.0; `"ligand-pooled"` reproduces the old behaviour. The general
lesson is worth more than the fix: **before trusting any pooled statistic, check how skewed the pool
is** -- `Counter(store._panel[cls].alleles).most_common(3)` is the whole check, and no test would
ever have caught this because every value was individually correct.

**A cache key of pure data cannot see a code change.** `predict.SCORER_EPOCH` is a hand-moved int in
the calibration-cache fingerprint, because 1.3.0 changed two scoring heads inside one released
version and the on-disk backgrounds cached before the change kept being served after it — which
silently moved a pinned test by 3.7x before it was caught. Bump it whenever the scoring code changes
what a head returns, and use a fresh `MHCMATCH_CALIBRATION_CACHE` for any run meant to *establish* a
number rather than iterate. `bench/results/calibration_cache_stale.md` has the measurement.

This file captures only *how we work in the repo*.

## Git worktrees — one worktree + branch per task

**We work in git worktrees.** Create a worktree (with its own branch) per feature/analysis rather than
switching branches in the main checkout, so parallel work never collides on `master`:

```zsh
git worktree add .claude/worktrees/<name> -b <name>     # isolated checkout on branch <name>
# ... work + commit inside .claude/worktrees/<name> ...
git worktree remove .claude/worktrees/<name>            # when the branch is consolidated
```

- One task → one worktree → one branch. Keep `<name>` short, kebab-case; reuse it as both dir and branch.
- `.claude/worktrees/` lives inside the repo dir but is **gitignored** (never commit nested checkouts).
- Consolidate finished results back and remove the worktree; merging `<name>` into `master` is a
  separate, deliberate step — `master` is never modified while parallel work is in flight.

## Four repos, four roles — never mix them

| local path | remote | holds | never holds |
|---|---|---|---|
| `~/vcs/code/mhcmatch` — **this one** | `antigenomics/mhcmatch`, **public** | library source, unit tests, sphinx docs, marimo notebooks, a handful of example images | benchmark harnesses, result tables, head-to-head comparisons, manuscript prose, publication figures |
| `~/vcs/projects/2026-mhcmatch-benchmark` | `repseq/2026-mhcmatch-code`, private → **released to reviewers** | every analysis script (`bench/`), every result table (`bench/results/*.md`), figure **generators**, `SOURCES.md` | library code, manuscript prose, any dataset |
| `~/vcs/manuscripts/2026-mhcmatch` | `repseq/2026-mhcmatch-ms`, private | manuscript, the **theory appendix** (`appendix/mhcmatch.tex`), every publication figure and generated LaTeX table | code that computes anything |
| `~/vcs/projects/2026-gamaleya-cancer` | `repseq/2026-gamaleya-cancer`, **private, stays private** | every run on real donors, and the donor key | anything that is ever copied outward under a name |

**The code repo is a reviewer artifact.** A reviewer clones `2026-mhcmatch-code`, runs
`mhcmatch bootstrap`, and every number regenerates. So what may be tracked there is exactly two
things: a **small metadata table** a script cannot derive, and a **result table**. Every dataset is
gitignored and bootstrapped from `~/hf/pmhc_data`; every feature frame is computed on the fly through
the `mhcmatch` CLI. A 28 MB cached parquet is neither of the two things, and it was tracked until
2026-08-23.

**The donor carve-out.** Runs on Gamaleya donors are how we know the release still works end to end
— de novo epitope counts, score distributions, rank movement between versions, human and mouse,
through Nextflow. Those runs happen in `2026-gamaleya-cancer` and stay there. What may cross into a
public or reviewer-facing repo is a **count or a distribution keyed by a donor code** (`D1`..`D8`,
`analysis/paired_tumor_norm/out/d24_donor_key.tsv` is the only place the codes meet the names) —
never a surname, never a peptide, never a genotype. De-identify at the moment you write the row.

Numbers flow one way: **benchmark → manuscript.** The manuscript never recomputes a number; it cites
one the benchmark repo produced and recorded. `bench/results/...` paths anywhere in this repo resolve
in the benchmark repo, not here.

**Everything in *this* repo bootstraps its data from HuggingFace.** Docs, notebooks and examples fetch
through `mhcmatch bootstrap` (the CLI command) / `Store` from `isalgo/pmhc_data` and friends — never from a hard-coded
`~/hf/...` or `~/vcs/projects/...` path, so a `pip install mhcmatch` user can run every example. Local
mirrors are for the benchmark repo, where the data is large and the analysis is one-off.

A test that asserts a *benchmark* number (an AUC, a head-to-head win) is research and belongs in the
benchmark repo. A test that asserts an *API contract* (a probability is in [0,1], a refactor is
score-identical) is a unit test and belongs here.

## This repo is public — nothing non-public goes in it

**mhcmatch is a public repo. Never commit non-public data to it — patient or sample data, private
cohort outputs, anything under NDA or consent.** Not in code, tests, fixtures, docs, notebooks,
commit messages, or a benchmark table. A git history is public forever once pushed, and `git rm`
does not remove a blob that has already shipped.

This is a property of the *repo*, not of any one path. There used to be a `bench/results/private/`
ignore rule standing in for the policy; it went stale the moment `bench/` moved out, and a path rule
was never the right shape anyway — it protects one directory and says nothing about the next one
someone invents.

- Private analyses belong in a private repo (e.g. the per-patient concordance runs, which live with
  the benchmark harness). Reference their *conclusions* here, not their contents.
- Public-by-construction data is fine: the `isalgo/pmhc_data` compendium, IEDB, IMGT, published
  benchmarks.
- Sample identifiers count as data. A surname plus an HLA genotype is identifying, and so is a
  surname on its own in a table of clinical-cohort measurements. **The way this rule actually gets
  broken is a per-donor results table pasted out of the private repo into a public `.md` to make a
  point** -- that is how eight surnames reached `ROADMAP.md` and a public `master` on 2026-08-21.
  De-identify at the moment you write the row, not at review: the counts carry the argument and the
  names never do.
- If you are unsure whether something is public, it is not. Ask.

The mirror-image rule, for the same reason: **shipped package data is never ignorable.**
`.gitignore` carries a `*.fasta` glob for fetched proteomes, and `src/mhcmatch/data/mhci_pseudo.fa`
escapes it only by its extension. A `!src/mhcmatch/data/**` negation is what actually holds — without
it, a model file named `.fasta` would be silently untracked and every fresh clone would ship a
package that cannot load its own models.

## Benchmarks — record the result, and scrutinise asymmetrically *on purpose*

**Every benchmark run that completes gets recorded. Never delete a result because we won it.**

- **Where we win:** report it as measured. Do not go hunting for a reason it might not count, and do
  not suppress it pending one. If something about the run is notable (an opponent below chance, a
  stratum of n=1), write it down *next to* the number as an observation and move on. Deciding what
  survives peer review is the author's call, not the run's.
- **Where we lose:** that is where the digging goes. Find the mechanism, fix it, re-run.

This is deliberate anti-symmetry. Scrutiny costs effort, so spending it only on the wins is how a
method gets talked down: every win acquires a caveat and every loss is taken at face value, and the
reported method is strictly worse than the real one. The bias to correct is ours, not the data's.

Precedent: `bench/results/compare_mhc2_mouse_hard_ligandbg.md` (a nine-cell sweep) was once deleted
for having an opponent score below chance. It is restored, with that fact recorded beside it.

## Data sources — one corpus, tuned per task

**Train on the whole corpus; do not filter it to make a benchmark look clean.** The general model is
fit on everything (EL, BA, in-silico) and beats broadly; per-task behaviour comes from **parameters**
(`background`, `footprint`, `register`, `h`, `tau`), not from a smaller training set. Binding-assay
peptides do bind — they are valid motif evidence.

Filtering is for **evaluation strata** (what a given number is *about*), never for training. That is
what `run_compare.py --el-only` is: it chooses which pairs may be positives, and the model behind it
is still the shipped one.

Provenance as a *model* parameter — an adjusted general model per source — was tested and is **not
needed**: the corpus-learned core-offset prior beats a uniform one by +0.010 on eluted-ligand queries
and +0.001 on binding-assay ones, i.e. it helps where boundaries carry information and is harmless
where they do not. Re-test if provenance ever reaches the pmhc schema; do not build the plumbing on
spec. See `bench/results/compare_mhc2_human_hard_ligandbg_elonly.md`.

## Reproducibility — everything ships from `mhcmatch build`

**Every artifact mhcmatch ships is rebuilt from the CLI, and the whole rebuild costs minutes. So
there is never a reason to run a stale one — always regenerate to the latest version.**

```zsh
mhcmatch build --check      # what is stale? builds nothing, exits 1 if any. This is what CI runs
mhcmatch build              # rebuild everything buildable in-process
mhcmatch build corpus -v    # one target, with per-step wall clock
```

- Builders live in `src/mhcmatch/_build.py`, **not** in `tools/`. A builder that only exists in a
  source checkout is one a wheel user cannot run, and one that gets forgotten. The
  `tools/build_anchor_models.py` / `tools/build_corpus_tables.py` shims are **gone** — a second name
  for one command is a second thing to keep current, and `publish.yml` had been calling the shim for
  one family while the other went unrebuilt. `tools/build_recognition.py` stays: it needs ESM2 and
  genuinely cannot run in-process, which is why `build` prints its command instead.
- **A new shipped artifact is not done until it has a `build` target and a `PROVENANCE.md` entry.**
- Inputs bootstrap from HuggingFace (`isalgo/pmhc_data`), never a local `~/hf/...` or
  `~/vcs/projects/...` path.
- **`aggregate_mhc1.json` (the EPIC scorer) has a real generator, and it still does not ship itself.**
  `bench/epic/fit.py --physchem rose_af5 --presentation binder` writes a **candidate** to
  `bench/epic/aggregate_mhc1.json` and deliberately does *not* copy it over
  `src/mhcmatch/data/aggregate_mhc1.json` -- what ships is the author's call, not the run's. `bench/run_epic.sh` is the whole chain that leads to the
  candidate -- deposits, features, Kd, Chowell, kernel grid, fit, corpus ladder, chemistry arms,
  selection, then the manuscript tables. Stage 0 gates on `mhcmatch --version` matching the
  checkout's `pyproject.toml`; that flag did not exist until 0.27.0, so the chain had never run
  past stage 0.

  **Because the copy is manual, the shipped artifact drifts, and `build --check` cannot see it** --
  it compares version stamps, and a hand-copied older fit stamped with the current version is
  current by that test. Measured 2026-08-23: the shipped file was the 06:15 fit against a frame
  whose upstream parquets were stale, and the chain at 08:50 produced BIC 4172.4 -> 4168.6, LOO
  mean 0.6602 -> 0.6654. After any run of the chain, diff the candidate against the shipped file
  before believing they agree.
  The binder pass is **not** the cost to optimise: `binder_score` on a warm allele is 82,201 pair/s,
  i.e. 4.4 s for all 363,324 pairs. The 1,314 s is per-allele calibrator construction -- 0.95 s cold
  x 203 alleles, per worker -- plus `Store.from_pmhc` at 4.5 s per worker per host. Batch the
  calibrator builds, not the peptide loop.
- On a version bump, regenerate rather than reasoning about whether it matters, then *verify* the
  scores did not move. Both existing builders are instrumented for this: `corpus_tables` prints
  `** MOVED **` for any cell that changed, and the anchor rebuild is checked against the previous
  file. Bump-only rebuilds have measured max |new − old| = 0.
- **Two version vocabularies, told apart by the shape of the value, not by the file extension.**
  A **model** version is an `int` — EPIC is `10`, the recognition heads are `2`, mimicry is `1`; a
  **package** version is dotted (`0.26.0`). `--check` compares the dotted ones to `__version__` and
  presence-checks the rest. Comparing a model version to a package version is a category error that
  reports every head stale at every release, which is why `.json` was once blanket-exempted — and
  the exemption cost the check *all* of them. `mimicry_fit.py` wrote `"0.12.0"` as a model version,
  the one file that made the shapes ambiguous; it now writes `1`, and the rule holds by construction.
- **`--check` covers all 27 shipped artifacts, not 11.** Until 2026-08-23 the `TARGETS` table listed
  only the anchors, the corpus tables and the recognition heads, so sixteen files — including
  `aggregate_mhc1.json` (EPIC itself) and `affinity_potts_mhc{1,2}.npz` (the source of `occupancy`)
  — could go missing or ship half-copied and nothing would say so. Entries whose generator lives in
  the benchmark repo carry a `None` builder and an `EXTERNAL` command taken from the artifact's own
  `generator` field or from `PROVENANCE.md`; `build` prints it. A target with **no** generator on
  record says so rather than printing an invented command.

**Why this is a rule.** 0.25.0 → 0.26.0 shipped three vendored `AnchorModel`s still stamped 0.25.0
*and* a `corpus_tables.npz` stamped 0.25.0. Only the models had a guard test, and it passed locally
because the editable install's metadata was also stale, so it compared the wrong value against
itself. CI caught one of the two; nothing would have caught the other.

## The scored columns are versioned, and the library carries both vocabularies

`aggregate_score` reads **only** the names in the artifact's `features` list, so `rank._finish`
supplies the union of every name a shipped or candidate artifact could ask for. That is what lets
one library score a v3 and a v4 artifact with no branch, and it is why adding a term is additive
rather than a migration:

| v3 name | v4 name | note |
|---|---|---|
| `binder` | `pres` | presentation head alone; `binder` folds the affinity rank in a second time |
| `occupancy` | `occupancy` | unchanged; `d_occupancy` is emitted and **not** fitted |
| -- | `d_occupancy`, `wt_absent` | emitted, measured, not fitted -- neither earned its parameter |
| `C_phys_rose` | `C_phys_buried` | same Rose scale, both keys computed |

**Do not delete an old name when a new one lands.** A recorded result cites the model version it was
produced under, and a registry that drops the old name cannot say what those numbers were.

**The corpus face and kernel are parameters, not constants.** `mimicry.face_kmers(mask=)` and
`contract(kernel=)` take `"slice"`/`"wildcard"` and any 20x20 or 21x21 array; `_build.SHIPPED_CORPUS`
names the `(k, mask)` a release commits to. A table is a pure function of
`(cls, comp, species, k, mask)` and the vendored key encodes all five -- so a wildcard-masked query
cannot silently index a sliced table. Hamming is kept only so pre-0.27 results reproduce; **the
corpus channels are BLOSUM62 from v4 on.**

## Class II — the levers, and the one that is not a lever

Four class-II knobs are measured against the shipped `register_em=2, tau=10, K=3, footprint=adaptive`
on both human arms (`bench/results/register_em_convergence_dp.md`,
`mhc2_register_frequency_gate.md`). **None dominates the shipped default, and two of the four optima
are mutually exclusive by construction** -- do not re-derive this:

| config | screening rare AUPRC | screening frequent AUPRC | screening frequent PPV@P |
|---|--:|--:|--:|
| shipped `em=2`, `tau=10` | 0.648 | 0.625 | 0.579 |
| `converge-frequent`, `tau=10` | 0.641 | **0.667** | 0.619 |
| `converge-frequent`, `tau=auto` | 0.639 | **0.668** | **0.629** |
| `em=2`, `tau=auto` | **0.689** | 0.616 | 0.575 |

Rare wants `em=2`; frequent wants `converge`. They do not compose, and the reason is structural: a
rare allele's motif is **67-77% borrowed** from its groove neighbours, and the neighbours are exactly
the alleles convergence moves. Gating the *borrower* -- its register, its mixture, its null, its
donor table, its `tau` -- was built (`register_em="converge-frequent"`) and recovers only about half
the rare screening loss. Getting both optima needs **routing by allele frequency at the model level**
(two fitted models), which is a product decision, not a parameter.

**`footprint=anchor` is the one setting that is never right on the predict path.** `MHC2_ANCHORS =
(1,4,6,9)` reaches only `Store._anchor_model` (`restriction`, `vote`); `build_scorer` ships
`adaptive`, which maps to all nine core positions. A benchmark arm left at `anchor` understates
mhcmatch -- `compare_mhc2_mouse_hard_ligandbg.md` did, on all nine cells. Reached from the other
side, `families=` (per-component gap placement) can only *subtract* from a 9-mer core and loses
monotonically in how many positions it masks.

## MixMHC2pred is the third rival, and it is the informative one

Installed at `~/work/academy/software/MixMHC2pred-2.1` (v2.1-beta1); adapter `bench/compare/mixmhc2.py`,
analysis `bench/mixmhc2/discordance_mix.py`, provenance `bench/mixmhc2/SOURCES.md`. It is worth more
than a second AUC because it is architecturally the same object we fit -- per-allele PWM mixtures with
an explicit binding core -- and it **reports which component fired**. Two facts that came only from it:

- **`SubSpec = -1` is reverse-orientation binding**, a mode `AnchorModel` cannot represent at all, and
  it is enriched in our misses. `AnchorModel(reverse=p)` marginalises the C-to-N reading with prior
  mass `p`; **ships off**, and `reverse=0.0` is bit-identical.
- **Allele names drop silently.** An unresolved allele produces no `%Rank_<allele>` column and the
  binary still exits 0, so `mix_allele` checks the shipped `PWMdef/` rather than trusting its own
  conversion. mhcmatch covers 47 of 47 human class-II panel alleles; MixMHC2pred ships PWMs for 42.

## Git flow & commits

- Branch flow: **feature → `dev` → `master`** (`ROADMAP.md` §7).
- End commit messages with the `Co-Authored-By` trailer. No PyPI release without explicit sign-off.
- Never fabricate citations — verify every DOI via a tool before adding it to `../../manuscripts/2026-mhcmatch/appendix/refs.bib`.

## Environment

- Repo-local `.venv` for the library. `environment.yml` is the heavier **notebooks** env
  (`mhcmatch-bench`: mmseqs2, graphviz, gnuplot, editable `../tcren-ms`); the `bench/` head-to-heads
  it was originally written for live in the benchmark repo now, along with their datasets at
  `~/hf/pmhc_data`.
