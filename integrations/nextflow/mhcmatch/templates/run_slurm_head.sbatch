#!/bin/bash -l
#SBATCH --job-name=mhcmatch-head
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=mhcmatch-head-%j.out
#
#   sbatch --partition <your-partition> run_slurm_head.sbatch
#
# The COHORT form. Nextflow's own process is what sbatch runs; it then submits one job per task, so
# samples run in parallel across the cluster instead of inside one allocation. Give the head job a
# small allocation and a long wall clock, because it mostly waits.
set -euo pipefail

# ---------------------------------------------------------------- EDIT THESE
ENV=mhcmatch
MODULE=/path/to/mhcmatch/integrations/nextflow/mhcmatch
REF=/shared/ref/mhcmatch
IN=/path/to/donor_files
OUT=$PWD/results
QUEUE=                             # REQUIRED. `--mhcmatch_slurm_queue` falls back to `normal`,
                                   # which is a common partition name and not a universal one, so
                                   # the default gets every task rejected at submit on a cluster
                                   # that does not have it. Check `sinfo`, and check its TIME LIMIT:
                                   # MHCMATCH_CASSETTE asks for 8 h under --screen and a 2 h queue
                                   # cannot give it
N0=8                               # per-allotype capacity. NOT read on this path -- both arms
                                   # run `cassette order`, which does not size, and
                                   # `cassette select` sizes by fixed -k. Only
                                   # `cassette build` reads it.
K=20                               # EPITOPES selected -- the construct-size lever here. The
                                   # cassette carries no more units than this, usually fewer.
TUMOR=SKCM
# -----------------------------------------------------------------------------
[ -n "$QUEUE" ] || { echo "set QUEUE: run 'sinfo -o \"%20P %5a %10l %6D %6c %10m\"' first" >&2; exit 2; }

# **`MHCMATCH_PMHC_DIR` is a READ override, not a download destination.** It is consulted first
# and used when the file is already there; when it is not, the fetch goes through
# `hf_hub_download`, which writes to the HuggingFace cache and ignores this variable entirely. So
# `HF_HOME` is the one that decides where ~250 MB of reference data actually lands -- without it
# that is your home directory, which on a shared cluster is the small quota this block exists to
# avoid. Measured on one cluster: 249 MB in $HOME with MHCMATCH_PMHC_DIR set and honoured.
export HF_HOME="$REF/hf"
export MHCMATCH_PMHC_DIR="$REF/pmhc_data"
export MHCMATCH_CALIBRATION_CACHE="$REF/calibration"
export TMPDIR="${TMPDIR:-$REF/tmp}"
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME="$REF/.nextflow"
mkdir -p "$HF_HOME" "$TMPDIR" "$NXF_HOME"

# **Each run gets its own launch directory.** Nextflow keeps its session cache under the directory
# it was launched from, so two runs started from the same one -- human and mouse side by side --
# collide on `Unable to acquire lock on session with ID ...` and the second dies at startup. It is
# also where `.nextflow.log` lands, which you want separate when something fails.
mkdir -p work_cohort && cd work_cohort

# **On THIS path the three exports above are not enough, and the three `--mhcmatch_*` flags below
# are what actually place the data.** `-profile slurm` loads slurm.config, whose `env {}` block
# *sets* the task environment from `params.mhcmatch_{pmhc_dir,calibration_cache,hf_home}` rather
# than inheriting this shell's. Each task is its own sbatch, so the exports reach only the head
# process. Pass `--mhcmatch_hf_home` and the tasks write to `$REF/hf`; omit it and every one of
# them falls back to the param default `${projectDir}/reference/hf` -- ~250 MB inside the module
# checkout, which is the one directory that may be read-only or on a home quota. The local-executor
# templates (run_human/run_mouse.sbatch) need no flags: they never load slurm.config, so tasks
# inherit the exported environment.

conda run -n "$ENV" --no-capture-output nextflow run "$MODULE/pipeline.nf" \
    -profile slurm \
    --indir   "$IN" \
    --outdir  "$OUT" \
    --mhcmatch_slurm_queue       "$QUEUE" \
    --mhcmatch_pmhc_dir          "$REF/pmhc_data" \
    --mhcmatch_calibration_cache "$REF/calibration" \
    --mhcmatch_hf_home           "$REF/hf" \
    --mhcmatch_vector_n0  $N0 \
    --mhcmatch_cassette_k $K \
    --mhcmatch_tumor      $TUMOR \
    -resume

# `-resume` is not optional in practice: MHCMATCH_CASSETTE --screen builds a whole-proteome index
# per register length, and a re-run without it repeats hours of work that has not changed.
echo "RUN OK"
