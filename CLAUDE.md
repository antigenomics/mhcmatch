# CLAUDE.md — working conventions for mhcmatch

**Authoritative context lives elsewhere — read it first:**
- [`ROADMAP.md`](ROADMAP.md) — the agent contract: what mhcmatch is, phase status, open loops.
- [`appendix/mhcmatch.tex`](appendix/mhcmatch.tex) — the method/statistics spec.

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

> **Benchmarks live in a separate repo.** `bench/` moved to
> [`2026-mhcmatch-benchmark`](https://github.com/antigenomics/2026-mhcmatch-benchmark) — the head-to-head harness, the `bench/results/*.md`
> tables referenced throughout, and their provenance notes. Paths like `bench/results/...`
> below resolve there, not here.

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
- Never fabricate citations — verify every DOI via a tool before adding it to `appendix/refs.bib`.

## Environment

- Repo-local `.venv` for the library; the `bench/` head-to-heads use the conda `mhcmatch-bench` env
  (mmseqs2, gnuplot, editable `../tcren-ms`) — see `environment.yml`. Datasets at `~/hf/pmhc_data`.
