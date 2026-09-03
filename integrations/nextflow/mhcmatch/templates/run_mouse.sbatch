#!/bin/bash -l
#SBATCH --job-name=mhcmatch-mouse
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
# **Sized from measurement, not headroom-on-headroom.** Measured on Aldan-3 2026-09-03, on the
# published 1.7.2, both arms with the screen ON: two human donors plus two mouse lines ran in
# **26:48** at **16.3 GB** peak RSS, and two mouse lines alone in **12:30** at **12.7 GB**. The
# longest single task was the whole-proteome screen at **701 s**; nothing else exceeded 186 s.
# These asked 48 GB / 8 h -- 2.9x the memory and 18x the wall clock actually used, which on a
# shared cluster is time spent queueing behind a quota for capacity the job never touches. Raise
# them for a genuinely larger cohort, and raise the TIME first: memory tracks how many tasks the
# local executor runs at once, not how many samples you have.
#SBATCH --output=mhcmatch-mouse-%j.out
#
#   sbatch --partition <your-partition> run_mouse.sbatch
#
# Differs from run_human.sbatch in exactly three places, and each one is a decision:
#
#   --genome GRCm39            the only thing that selects `--species mouse`. There is no separate
#                              parameter, so mhcmatch follows the assembly the rest of a pipeline
#                              already runs on
#   --alleles / --alleles_mhc2 a literal panel instead of a typing file. Not a shortcut: an inbred
#                              line's H-2 haplotype is a property of the LINE, so there is nothing
#                              to type and nothing to infer
#   no --mhcmatch_tumor        the expression contexts are TCGA study codes and there is no mouse
#                              equivalent. Setting one scores mouse candidates against a human
#                              transcriptome's abundance floor
set -euo pipefail

# ---------------------------------------------------------------- EDIT THESE
ENV=mhcmatch
MODULE=/path/to/mhcmatch/integrations/nextflow/mhcmatch
REF=/shared/ref/mhcmatch
IN=/path/to/mouse_files
OUT=$PWD/results_mouse
N0=8                               # per-allotype capacity. NOT read on this path -- both arms
                                   # run `cassette order`, which does not size, and
                                   # `cassette select` sizes by fixed -k. Only
                                   # `cassette build` reads it.
K=20                               # EPITOPES selected -- the construct-size lever here. The
                                   # cassette carries no more units than this, usually fewer.
# BALB/c (H-2d) below. C57BL/6 (H-2b) is 'H2-K*b,H2-D*b' and 'H-2-IAb'. All three spellings resolve
# -- `H2-K*d`, `H-2Kb`, `I-Ab` -- so pass whatever your own tables carry.
ALLELES_MHC1='H2-K*d,H2-D*d,H2-L*d'
ALLELES_MHC2='H-2-IAd,H-2-IEd'
                      # ^ passed TWICE below: once as the class-II panel to score against, and once
                      # as the recipient's allotypes for the cassette map. Without the second the map
                      # is class I only and `self_help` -- does this unit's CD8 epitope have CD4 help
                      # from the SAME unit -- is not computed
# P(an allotype is live) in the response model behind the quota. 0.999 for mouse against 0.95 for
# human is what the shipped bundles used: a stated design parameter, not a fitted one. Measure your
# own with mhcmatch.portfolio.betabinom_rho before trusting either.
BLOCK_LIVE=0.999   # P(allotype alive) behind a QUOTA. Inert unless you also set
                                   # --mhcmatch_vector_quota: MHCMATCH_CASSETTE emits
                                   # --block-live only alongside --quota, and drops it
                                   # otherwise. `cassette select` uses the separate
                                   # --mhcmatch_cassette_block_live (default 1.0).
# -----------------------------------------------------------------------------

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
mkdir -p work_mouse && cd work_mouse

conda run -n "$ENV" --no-capture-output python -c "import mhcmatch; print('mhcmatch', mhcmatch.__version__)"

conda run -n "$ENV" --no-capture-output nextflow run "$MODULE/pipeline.nf" \
    --indir   "$IN" \
    --outdir  "$OUT" \
    --mode    both \
    --genome  GRCm39 \
    --alleles      "$ALLELES_MHC1" \
    --alleles_mhc2 "$ALLELES_MHC2" \
    --mhcmatch_vector_map_alleles_mhc2 "$ALLELES_MHC2" \
    --mhcmatch_vector_n0         $N0 \
    --mhcmatch_cassette_k        $K \
    --mhcmatch_vector_block_live $BLOCK_LIVE \
    -resume

echo "RUN OK"
