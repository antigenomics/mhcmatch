# mhcmatch epitope-prediction image, for the -profile docker path of the module beside it.
#
# Packages `mhcmatch predict` + its seqtree C++ core + the reference ligand panel, so a process can
# run:  mhcmatch predict <fasta> --alleles <...> --cls <mhc1|mhc2> --scored-csv <out> --native <out>
# with no host data mounts. The panel is bootstrapped from the public HF dataset isalgo/pmhc_data at
# BUILD time (baked into the image's huggingface_hub cache), so runtime needs no network -- which
# matters on offline compute nodes.
#
# Build (no data staging needed -- the panel is auto-fetched):
#   docker build -t <registry>/mhcmatch:0.4.1 integrations/nextflow/mhcmatch/
#
# mhcmatch and seqtree come from PyPI (pinned via MHCMATCH_VERSION). The build tools are a fallback
# for platforms without a prebuilt seqtree wheel; on ARM/x86_64 linux the wheel is used directly.

FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build \
    && rm -rf /var/lib/apt/lists/*

ARG MHCMATCH_VERSION=0.4.1
RUN pip install --no-cache-dir "mhcmatch==${MHCMATCH_VERSION}"

# Bootstrap the reference ligand panel from the public HF dataset (both tiers, ~12 MB + ~4 MB) into
# the image's huggingface_hub cache, so from_pmhc() resolves it offline at runtime.
RUN mhcmatch bootstrap

# Sanity: the CLI and panel resolve at build time (fails the build early if the panel is wrong).
RUN mhcmatch --help >/dev/null && \
    python -c "from mhcmatch import Store; Store.from_pmhc(tier='shortlist', species='human', classes=('mhc1',))"

ENTRYPOINT ["mhcmatch"]
CMD ["--help"]
