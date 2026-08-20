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
- Sample identifiers count as data. A surname plus an HLA genotype is identifying.
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

## Git flow & commits

- Branch flow: **feature → `dev` → `master`** (`ROADMAP.md` §7).
- End commit messages with the `Co-Authored-By` trailer. No PyPI release without explicit sign-off.
- Never fabricate citations — verify every DOI via a tool before adding it to `../../manuscripts/2026-mhcmatch/appendix/refs.bib`.

## Environment

- Repo-local `.venv` for the library. `environment.yml` is the heavier **notebooks** env
  (`mhcmatch-bench`: mmseqs2, graphviz, gnuplot, editable `../tcren-ms`); the `bench/` head-to-heads
  it was originally written for live in the benchmark repo now, along with their datasets at
  `~/hf/pmhc_data`.
