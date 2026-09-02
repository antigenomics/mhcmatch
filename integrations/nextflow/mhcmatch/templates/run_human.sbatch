#!/bin/bash -l
#SBATCH --job-name=mhcmatch
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mhcmatch-%j.out
#
#   sbatch --partition <your-partition> run_human.sbatch
#
# The LOCAL executor inside one allocation, which is what a handful of samples wants. For a real
# cohort use run_slurm_head.sbatch: it submits one job per task, which is the point of a cluster.
#
# Run setup.sbatch first. Nothing here installs anything -- two runs sharing one conda env race on
# `pip install --force-reinstall` and one clobbers the other's entry point mid-write.
set -euo pipefail

# ---------------------------------------------------------------- EDIT THESE
ENV=mhcmatch
MODULE=/path/to/mhcmatch/integrations/nextflow/mhcmatch    # this repo, checked out
REF=/shared/ref/mhcmatch
IN=/path/to/donor_files            # the naming contract is in ../README.md
OUT=$PWD/results
TUMOR=SKCM                         # TCGA study code. `mhcmatch expression --list-contexts` prints
                                   # the 19 TCGA<->GTEx pairings. Leaving it unset answers "is this
                                   # gene expressed anywhere" instead of "in this tumour" -- set it
N0=8                               # per-allotype capacity. REQUIRED and with no default on
                                   # purpose: nothing in the public record fits it, so the value is
                                   # yours to defend and it is recorded in the output
K=20                               # units manufactured. A different question from N0
# -----------------------------------------------------------------------------

export MHCMATCH_PMHC_DIR="$REF/pmhc_data"
export MHCMATCH_CALIBRATION_CACHE="$REF/calibration"
export TMPDIR="${TMPDIR:-$REF/tmp}"
export NXF_ANSI_LOG=false
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME="$REF/.nextflow"
mkdir -p "$TMPDIR" "$NXF_HOME" logs

conda run -n "$ENV" --no-capture-output python -c "import mhcmatch; print('mhcmatch', mhcmatch.__version__)"

conda run -n "$ENV" --no-capture-output nextflow run "$MODULE/pipeline.nf" \
    --indir   "$IN" \
    --outdir  "$OUT" \
    --mode    both \
    --mhcmatch_vector_n0  $N0 \
    --mhcmatch_cassette_k $K \
    --mhcmatch_tumor      $TUMOR \
    -resume

# `--mhcmatch_vector_screen` defaults to TRUE and should stay that way. With it off, NO safety
# check runs at all and the cassette carries whatever it was handed. It is what the 48 GB and the
# 8 h above are for: one whole-proteome index per register length, ~12 GB peak each. Adding
# `--mhcmatch_vector_screen false` makes the run minutes instead of hours and is for iterating on
# the plumbing -- never for a cassette anyone will manufacture.
echo "RUN OK"
