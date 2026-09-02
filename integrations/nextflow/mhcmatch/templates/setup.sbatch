#!/bin/bash -l
#SBATCH --job-name=mhcmatch-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/mhcmatch-setup-%j.out
#
# Run this ONCE, before the first pipeline run, and never again.
#
#   sbatch --partition <your-partition> setup.sbatch
#
# It is a JOB and not a login-node command on purpose: `bootstrap --reference` downloads and
# unpacks ~115 MB, and login nodes are shared, often two-core, and on some clusters a guard kills
# anything that looks like compute there. Compute nodes usually have outbound HTTPS.
#
# `#!/bin/bash -l` and `conda run --no-capture-output`, never `conda activate`: conda's shell hook
# is only loaded by a login shell, and `activate` in a non-interactive batch script silently leaves
# you on the system interpreter.
set -euo pipefail

# ---------------------------------------------------------------- EDIT THESE
ENV=mhcmatch          # conda env name
REF=/shared/ref/mhcmatch
                      # the shared reference + calibration root. A filesystem EVERY compute node can
                      # see, and a large one -- not a home quota, which on a shared cluster is
                      # typically small, untracked and never cleaned up
TYPING=""             # optional but worth it: one typing file. Prints how many alleles resolved,
                      # which is the five-second check that a run is not about to score against an
                      # empty panel
WHEELHOUSE=""         # optional: a directory of .whl files, for a node whose egress cannot reach
                      # PyPI. See the note further down
# -----------------------------------------------------------------------------

export MHCMATCH_PMHC_DIR="$REF/pmhc_data"
export MHCMATCH_CALIBRATION_CACHE="$REF/calibration"
export TMPDIR="${TMPDIR:-$REF/tmp}"
mkdir -p "$MHCMATCH_PMHC_DIR" "$MHCMATCH_CALIBRATION_CACHE" "$TMPDIR" logs

# **This is the only place that installs.** Two runs sharing one conda env race on
# `pip install --force-reinstall` and one clobbers the other's entry point mid-write.
if [ -n "$WHEELHOUSE" ]; then
    # A compute node's egress may not reach PyPI even when it reaches HuggingFace -- measured on one
    # cluster as a read-timeout against pypi.org after four retries, while `bootstrap` succeeded from
    # the same node. Build the wheelhouse where the network works:
    #   pip download --dest wheels --platform manylinux2014_x86_64 --python-version 3.12 \
    #       --implementation cp --only-binary=:all: mhcmatch
    conda run -n "$ENV" --no-capture-output pip install -q --no-index --find-links "$WHEELHOUSE" --force-reinstall mhcmatch
else
    conda run -n "$ENV" --no-capture-output pip install -q --upgrade mhcmatch
fi

echo "== versions =="
conda run -n "$ENV" --no-capture-output mhcmatch --version
conda run -n "$ENV" --no-capture-output nextflow -version 2>&1 | grep -i version || true

echo "== references (~115 MB: ligand panel, thymic/viral sets, neoantigen deposits, expression) =="
conda run -n "$ENV" --no-capture-output mhcmatch bootstrap --reference
conda run -n "$ENV" --no-capture-output mhcmatch bootstrap --proteome human,mouse

# A five-second check that is worth more than it looks: every HLA typer writes the G-group form
# (`A*01:01:01G`), and an allele that resolves to nothing is dropped SILENTLY one layer down --
# so a run against an empty panel exits 0. This says how many resolved, and names every drop.
if [ -n "$TYPING" ]; then
    echo "== allele resolution on $TYPING =="
    conda run -n "$ENV" --no-capture-output mhcmatch alleles "$TYPING" --cls mhc1
    conda run -n "$ENV" --no-capture-output mhcmatch alleles "$TYPING" --cls mhc2
fi

echo "SETUP OK"
