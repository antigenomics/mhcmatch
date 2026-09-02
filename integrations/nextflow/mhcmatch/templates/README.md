# Templates — copy these, edit the block at the top, submit

Four SLURM scripts. Each has one clearly-marked `EDIT THESE` block near the top and nothing
cluster-specific below it. Everything they need is in this repository; nothing else has to be
fetched.

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

**`--mhcmatch_vector_n0` has no default on purpose.** Per-allotype capacity is not fitted by
anything in the public record, so the value is yours to defend — and it is recorded in the output so
a reader can see which one you chose. It is a different question from `--mhcmatch_cassette_k`:
`k` is how many units are *manufactured*, `n0` is how many the recipient's allotypes can *carry*.

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
