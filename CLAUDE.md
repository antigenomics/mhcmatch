# CLAUDE.md — working conventions for mhcmatch

**Authoritative context lives elsewhere — read it first:**
- [`ROADMAP.md`](ROADMAP.md) — the agent contract: what mhcmatch is, phase status, open loops.
- `../../manuscripts/2026-mhcmatch/appendix/mhcmatch.tex` — the method/statistics spec (manuscript repo).

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

## Three repos, three roles — never mix them

| repo | holds | never holds |
|---|---|---|
| `~/vcs/code/mhcmatch` — **this one** | library source, unit tests, sphinx docs, marimo notebooks, a handful of example images | benchmark harnesses, result tables, head-to-head comparisons, manuscript prose, publication figures |
| `~/vcs/projects/2026-mhcmatch-benchmark` | every analysis script (`bench/`), every result table (`bench/results/*.md`), figure **generators**, `SOURCES.md` | library code, manuscript prose |
| `~/vcs/manuscripts/2026-mhcmatch` | manuscript, the **theory appendix** (`appendix/mhcmatch.tex`), every publication figure and generated LaTeX table | code that computes anything |

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
  source checkout is one a wheel user cannot run, and one that gets forgotten. `tools/build_*.py`
  are shims kept only because `PROVENANCE.md` names those paths.
- **A new shipped artifact is not done until it has a `build` target and a `PROVENANCE.md` entry.**
- Inputs bootstrap from HuggingFace (`isalgo/pmhc_data`), never a local `~/hf/...` or
  `~/vcs/projects/...` path.
- **`aggregate_mhc1.json` (the EPIC scorer) has no `build` target yet, and this is the open loop.**
  Not for want of data: all four labelled deposits it is fitted on are already published at
  `isalgo/pmhc_data/neoantigens/` (`neoag_tested.tsv.gz` 321,825 rows, `neoantigens_tested_peptides.tsv.gz`
  31,804, `neoag_tested_mmu.tsv.gz` 866, `neoag_tested_hsa.tsv.gz` 414), and of the fit frame's 47
  columns only **28 are irreducible** — 4.7 MB of parquet — while the other 19 are computed by
  mhcmatch itself. What is missing is that the assembly, the feature build and the fit live in the
  benchmark repo (`corpus_grand.py` -> `grand_corpus.py` -> `grand_ship.py`) and the artifact is
  hand-copied across. That hand-copy is how the GRAND -> EPIC rename reached the artifact but not its
  generator. Closing it needs the chain moved into `_build.py` **and** the binder pass batched: it
  was recorded at 338,319 of 363,324 pairs in 1,314 s, which busts the budget until it gets the same
  one-batched-call treatment that took `store_binder` 15.8x.
- On a version bump, regenerate rather than reasoning about whether it matters, then *verify* the
  scores did not move. Both existing builders are instrumented for this: `corpus_tables` prints
  `** MOVED **` for any cell that changed, and the anchor rebuild is checked against the previous
  file. Bump-only rebuilds have measured max |new − old| = 0.
- Only the pickles and the `.npz` stamp `__version__`. A `.json` artifact's `version` is a **model**
  version — EPIC is `3`, the recognition heads are `2` — so `--check` validates those for presence
  only. Comparing the two is a category error that reports every head stale at every release.

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

## Git flow & commits

- Branch flow: **feature → `dev` → `master`** (`ROADMAP.md` §7).
- End commit messages with the `Co-Authored-By` trailer. No PyPI release without explicit sign-off.
- Never fabricate citations — verify every DOI via a tool before adding it to `../../manuscripts/2026-mhcmatch/appendix/refs.bib`.

## Environment

- Repo-local `.venv` for the library. `environment.yml` is the heavier **notebooks** env
  (`mhcmatch-bench`: mmseqs2, graphviz, gnuplot, editable `../tcren-ms`); the `bench/` head-to-heads
  it was originally written for live in the benchmark repo now, along with their datasets at
  `~/hf/pmhc_data`.
