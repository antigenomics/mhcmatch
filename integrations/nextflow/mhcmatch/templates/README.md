# Templates — copy these, edit the block at the top, submit

Four SLURM scripts. Each has one clearly-marked `EDIT THESE` block near the top and nothing
cluster-specific below it. Everything they need is in this repository; nothing else has to be
fetched.

**The only prerequisite is a working `conda`** (miniforge or miniconda). `setup.sbatch` creates
the environment itself — python 3.12 plus nextflow — installs the pinned mhcmatch release into it,
and refuses to continue if the version it got is not the one it asked for. You do not need to make
an env, and you do not need to install nextflow separately.

```bash
git clone https://github.com/antigenomics/mhcmatch.git
cp mhcmatch/integrations/nextflow/mhcmatch/templates/*.sbatch .
$EDITOR setup.sbatch && sbatch --partition <your-partition> setup.sbatch
```

Or take one file:

```bash
curl -O https://raw.githubusercontent.com/antigenomics/mhcmatch/master/integrations/nextflow/mhcmatch/templates/run_human.sbatch
```

| file | when | asks for |
|---|---|---|
| `setup.sbatch` | **once**, before anything else | 4 cpu · 16 GB · 1 h |
| `run_human.sbatch` | a handful of samples, local executor in one allocation | 8 cpu · 48 GB · 8 h |
| `run_mouse.sbatch` | the same, for an inbred line | 8 cpu · 48 GB · 8 h |
| `run_slurm_head.sbatch` | a **cohort** — one job per task, across the cluster | 2 cpu · 8 GB · 48 h |

### Every `EDIT THESE` variable

| variable | in | default | what it is |
|---|---|---|---|
| `ENV` | all four | `mhcmatch` | conda env name. `setup.sbatch` **creates** it if absent |
| `VERSION` | setup | `1.7.3` | the release to install. Pinned, and asserted after install |
| `REF` | all four | `/shared/ref/mhcmatch` | shared reference + calibration root. A filesystem every compute node can see |
| `TYPING` | setup | *(empty)* | optional: one typing file, to print how many alleles resolve |
| `WHEELHOUSE` | setup | *(empty)* | optional: a directory of `.whl`s, for a node with no PyPI egress |
| `MODULE` | the three run scripts | `/path/to/…` | **path to this checkout's `integrations/nextflow/mhcmatch`** |
| `IN` | the three run scripts | `/path/to/…` | the input directory; naming contract in `../README.md` |
| `OUT` | the three run scripts | `$PWD/results…` | where results are published |
| `TUMOR` | human, head | `SKCM` | TCGA study code for expression context. **Do not set it for mouse** |
| `K` | the three run scripts | `20` | **epitopes** selected — the construct-size lever on this path |
| `N0` | the three run scripts | `8` | per-allotype capacity. **Not read on this path** (both arms run `cassette order`); only `cassette build` uses it |
| `QUEUE` | head only | *(empty)* | **REQUIRED.** `--mhcmatch_slurm_queue` falls back to `normal`, which many clusters do not have |
| `ALLELES_MHC1` / `ALLELES_MHC2` | mouse only | H-2d | the inbred line's haplotype; a literal panel replaces the typing file |
| `BLOCK_LIVE` | mouse only | `0.999` | P(allotype alive) **behind a quota** — inert unless `--mhcmatch_vector_quota` is also set |

`../README.md` is the module contract — every process's inputs and outputs, every parameter, the
input file-naming rules. `docs/pipeline.rst` is the same thing as prose.

## The three that are easy to get wrong

**`setup.sbatch` runs once, and it is the only thing that installs.** Two runs sharing one conda
env race on `pip install --force-reinstall` and one clobbers the other's entry point mid-write. The
run scripts only assert that the import works.

**`--mhcmatch_slurm_queue` has no safe default.** It falls back to `normal`, a common partition
name and not a universal one, so on a cluster without it every task is rejected at submit. Run
`sinfo -o '%20P %5a %10l %6D %6c %10m'` and check the **time limit** too: `MHCMATCH_CASSETTE` asks
for 8 h under `--screen`, which a 2 h queue cannot give it.

**`--mhcmatch_vector_n0` is not read on the `pipeline.nf` path at all.** Both arms run
`cassette order`, which does not size, and `cassette select` sizes by fixed `-k` -- so
`--mhcmatch_cassette_k` is the construct-size commitment here. `--n0` is read by
`cassette build`, which is what `subworkflows/mhcmatch.nf` runs, and there it
**has no default on purpose.** Per-allotype capacity is not fitted by
anything in the public record, so the value is yours to defend — and it is recorded in the output so
a reader can see which one you chose. It is a different question from `--mhcmatch_cassette_k`:
`k` is how many **epitopes are selected** — the construct carries fewer *units* than that,
because several epitopes can share one 27-mer window and the screen withdraws some (20 → 15
→ 11 on one measured donor; see `-k` counts epitopes, not manufactured units in
`../README.md`). `n0` is how many the recipient's allotypes can *carry*.

**Two runs cannot share a launch directory.** Nextflow keeps its session cache under the directory
it was launched from, so `run_human.sbatch` and `run_mouse.sbatch` started side by side would
collide on `Unable to acquire lock on session with ID ...` and the second would die at startup.
Each script `cd`s into its own `work_*` subdirectory for that reason — which is also where its
`.nextflow.log` lands, and you want those separate when something fails.

## Two things worth reading before you change a default

**`--mhcmatch_vector_screen` is on, and off means off.** With it false, **no** safety check runs at
all and the cassette carries whatever it was handed. It is what the 48 GB and the 8 h are for. Turn
it off to iterate on plumbing; never for a cassette anyone will manufacture.

**A compute node's egress may not reach PyPI even when it reaches HuggingFace.** Measured on one
cluster as a read-timeout against `pypi.org` after four retries, while `mhcmatch bootstrap` from the
*same* node succeeded. `setup.sbatch` has a `WHEELHOUSE` knob for it:

```bash
# where the network works:
pip download --dest wheels --platform manylinux2014_x86_64 --python-version 3.12 \
    --implementation cp --only-binary=:all: mhcmatch
# then set WHEELHOUSE=/path/to/wheels in setup.sbatch
```

## Checking it worked

`setup.sbatch` ends in `SETUP OK` and the run scripts in `RUN OK`. Set `TYPING=/path/to/typing.tsv`
before submitting `setup.sbatch` and it also prints how many alleles resolved — worth the five
seconds, because an allele that resolves to nothing is dropped **silently** one layer down, and a
run against an empty panel exits 0:

```
# dropped 6 name(s) that resolve to no pseudosequence: E*01:01, E*01:03, F*01:01, G*01:01
# 6 mhc1 allele(s) from 26 typed name(s)
```

Six is a plausible class-I count; zero means the typing file is not being read the way you think.
