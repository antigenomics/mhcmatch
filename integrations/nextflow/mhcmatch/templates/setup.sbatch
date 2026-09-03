#!/bin/bash -l
#SBATCH --job-name=mhcmatch-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=mhcmatch-setup-%j.out
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
ENV=mhcmatch          # conda env name. Created below if it does not exist yet
VERSION=1.9.0         # the mhcmatch release to install. PIN IT -- the pipeline calls
                      # subcommands (`mhcmatch alleles`) and flags (`rank --passthrough`)
                      # that do not exist before 1.8.0, and an unpinned `--upgrade` on a
                      # cluster with a stale index resolves to whatever it already has
REF=/shared/ref/mhcmatch
                      # the shared reference + calibration root. A filesystem EVERY compute node can
                      # see, and a large one -- not a home quota, which on a shared cluster is
                      # typically small, untracked and never cleaned up
TYPING=""             # optional but worth it: one typing file. Prints how many alleles resolved,
                      # which is the five-second check that a run is not about to score against an
                      # empty panel
WHEELHOUSE=""         # optional: a directory of .whl files, for a node whose egress cannot reach
                      # PyPI. See the note further down
STAGE_INDEX=""        # optional, and empty is right unless you enable the safety screen or the
                      # mimicry annotation -- BOTH SHIP OFF, and neither reads this. If you do turn
                      # one on, set e.g. "human:8|9|10|11" (GB per length) so the tasks load a
                      # prebuilt index instead of every one of them building its own at once
# -----------------------------------------------------------------------------

# **This script is the whole prerequisite.** It used to assume `$ENV` already existed with
# nextflow in it, which is true on a cluster where someone had already done it by hand and false
# for the colleague this template is written for. Creating an env that exists is an error, so the
# existence check comes first.
if ! conda env list | awk '{print $1}' | grep -qx "$ENV"; then
    echo "== creating conda env '$ENV' (python 3.12 + nextflow) =="
    # `--override-channels` is not tidiness. Without it conda also consults Anaconda's own
    # `defaults` channels, which are behind a Terms-of-Service acceptance gate; on a machine where
    # that has not been accepted the solve dies with "Terms of Service have not been accepted" and
    # creates nothing. conda-forge and bioconda carry everything needed and have no such gate.
    conda create -y -n "$ENV" --override-channels \
        -c conda-forge -c bioconda python=3.12 nextflow
    # conda can report success having created nothing; check rather than trust.
    conda env list | awk '{print $1}' | grep -qx "$ENV" || {
        echo "ERROR: conda reported success but '$ENV' does not exist" >&2; exit 1; }
fi

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
mkdir -p "$HF_HOME" "$MHCMATCH_PMHC_DIR" "$MHCMATCH_CALIBRATION_CACHE" "$TMPDIR"

# **This is the only place that installs.** Two runs sharing one conda env race on
# `pip install --force-reinstall` and one clobbers the other's entry point mid-write.
if [ -n "$WHEELHOUSE" ]; then
    # A compute node's egress may not reach PyPI even when it reaches HuggingFace -- measured on one
    # cluster as a read-timeout against pypi.org after four retries, while `bootstrap` succeeded from
    # the same node. Build the wheelhouse where the network works:
    #   pip download --dest wheels --platform manylinux2014_x86_64 --python-version 3.12 \
    #       --implementation cp --only-binary=:all: mhcmatch
    conda run -n "$ENV" --no-capture-output pip install -q --no-index --find-links "$WHEELHOUSE" --force-reinstall "mhcmatch==$VERSION"
else
    conda run -n "$ENV" --no-capture-output pip install -q --upgrade "mhcmatch==$VERSION"
fi

echo "== versions =="
conda run -n "$ENV" --no-capture-output mhcmatch --version
# A wrong version does not fail loudly on its own: the pipeline's first process calls
# `mhcmatch alleles`, which simply does not exist before 1.8.0, and the run dies several minutes
# later inside a task log. Check it here, where the message can say what actually happened.
conda run -n "$ENV" --no-capture-output mhcmatch --version 2>&1 | grep -q "$VERSION" || {
    echo "ERROR: expected mhcmatch $VERSION; the pipeline needs >=1.8.0 for \`alleles\` and \`rank --passthrough\`" >&2
    exit 1
}
conda run -n "$ENV" --no-capture-output nextflow -version 2>&1 | grep -i version || true

echo "== references (~115 MB: ligand panel, thymic/viral sets, neoantigen deposits, expression) =="
conda run -n "$ENV" --no-capture-output mhcmatch bootstrap --reference
conda run -n "$ENV" --no-capture-output mhcmatch bootstrap --proteome human,mouse

# **Only if you turn on the safety screen or the mimicry annotation.** Both ship OFF, and with
# them off nothing below is read -- so the default setup skips this entirely and the pipeline still
# runs end to end in minutes.
#
# Turn either on and a whole-proteome index has to exist. It is a pure function of published data,
# and its in-memory cache is per PROCESS -- so without staging it here, every task in a fan-out
# misses at the same moment and builds its own: measured on Aldan-3 2026-09-03, ~194 s per mimicry
# task against 65.0 s warm, and it was the dominant stage of the whole run. Staged once into $REF,
# every node loads it in well under a second.
#
# Fetched from the published dataset where available and built locally otherwise; the command says
# which, per length, because the two differ by two orders of magnitude. GB per length, so ask only
# for what you need -- `8|9|10|11` is the class-I ladder.
if [ -n "$STAGE_INDEX" ]; then
    echo "== prebuilt indexes (GB-scale; only needed for --screen / --mimicry) =="
    conda run -n "$ENV" --no-capture-output mhcmatch bootstrap --index "$STAGE_INDEX"
fi

# A five-second check that is worth more than it looks: every HLA typer writes the G-group form
# (`A*01:01:01G`), and an allele that resolves to nothing is dropped SILENTLY one layer down --
# so a run against an empty panel exits 0. This says how many resolved, and names every drop.
if [ -n "$TYPING" ]; then
    echo "== allele resolution on $TYPING =="
    conda run -n "$ENV" --no-capture-output mhcmatch alleles "$TYPING" --cls mhc1
    conda run -n "$ENV" --no-capture-output mhcmatch alleles "$TYPING" --cls mhc2
fi

echo "SETUP OK"
