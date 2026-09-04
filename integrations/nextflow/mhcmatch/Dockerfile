# mhcmatch epitope-prediction image, for the -profile docker path of the module beside it.
#
# Packages the mhcmatch CLI + its seqtree C++ core + every reference the module's nine processes
# read, so predict / rank / neoag / mimicry / cassette / cassette score all run with no host data mounts. Everything is
# staged from the public HF dataset isalgo/pmhc_data at BUILD time, into the image's
# huggingface_hub cache, so runtime needs no network -- which matters on offline compute nodes.
#
# Build (no data staging needed -- the references are auto-fetched):
#   docker build -t <registry>/mhcmatch:1.11.0 integrations/nextflow/mhcmatch/
#
# mhcmatch and seqtree come from PyPI (pinned via MHCMATCH_VERSION). The build tools are a fallback
# for platforms without a prebuilt seqtree wheel; on ARM/x86_64 linux the wheel is used directly.

FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build \
    && rm -rf /var/lib/apt/lists/*

ARG MHCMATCH_VERSION=1.11.0
RUN pip install --no-cache-dir "mhcmatch==${MHCMATCH_VERSION}"

# Bootstrap the reference ligand panel from the public HF dataset (both tiers, ~12 MB + ~4 MB) into
# the image's huggingface_hub cache, so from_pmhc() resolves it offline at runtime.
#
# `--reference` additionally stages the known-epitope sets, the mimicry references and the expression
# tables (~115 MB), which `rank`, `neoag` and `mimicry` all read. It is not optional now that those
# processes exist: without it they reach for HuggingFace from a compute node and fail there, at the
# worst possible time, rather than here at build time.
RUN mhcmatch bootstrap --reference

# **The proteomes, which `--reference` does NOT stage and the safety screen cannot run without.**
# `cassette build --screen` calls `fetch_proteome(species)` to build its whole-proteome window
# index, and the screen is ON by default -- so without this line the one process that most needs to
# be offline reaches for HuggingFace from a compute node, hours into a run, which is the exact
# failure the paragraph above says this image exists to prevent. Both species, because the module
# routes `--genome GRCm39` to `--species mouse`. `templates/setup.sbatch` stages the same two.
RUN mhcmatch bootstrap --proteome human,mouse

# Sanity: the CLI, the panel AND the proteome resolve at build time, so a missing bootstrap fails
# the build rather than a run. The proteome check is `fetch_proteome`, the same call the screen
# makes -- checking only the panel is what let the missing `--proteome` above go unnoticed.
RUN mhcmatch --help >/dev/null && \
    python -c "from mhcmatch import Store; Store.from_pmhc(tier='shortlist', species='human', classes=('mhc1',))" && \
    python -c "from mhcmatch.store import fetch_proteome; [fetch_proteome(s) for s in ('human', 'mouse')]"

ENTRYPOINT ["mhcmatch"]
CMD ["--help"]
