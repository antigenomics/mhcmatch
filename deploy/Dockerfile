# mhcmatch epitope-prediction image for the Gamaleya nextflow_vaccine pipeline.
#
# Packages `mhcmatch predict` + its seqtree C++ core + the reference ligand panel, so a process can
# run:  mhcmatch predict <fasta> --alleles <...> --cls <mhc1|mhc2> --scored-csv <out> --native <out>
# with no host data mounts. The panel is bootstrapped from the public HF dataset isalgo/pmhc_data at
# BUILD time (baked into the image's huggingface_hub cache), so runtime needs no network — which
# matters on offline compute nodes.
#
# Build (no data staging needed — the panel is auto-fetched):
#   docker build -t <registry>/mhcmatch:<ver> -f deploy/Dockerfile deploy/
#
# NOTE for review: pins default to the public antigenomics GitHub refs; override with --build-arg to
# pin an exact commit/tag, or swap to an internal mirror. Building seqtree needs a C++ toolchain.

FROM python:3.12-slim AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
        git build-essential cmake ninja-build \
    && rm -rf /var/lib/apt/lists/*

# seqtree (C++ core + bindings) then mhcmatch. Pin these for reproducible images.
ARG SEQTREE_REF=master
ARG MHCMATCH_REF=master
RUN pip install --no-cache-dir "git+https://github.com/antigenomics/seqtree.git@${SEQTREE_REF}"
RUN pip install --no-cache-dir "git+https://github.com/antigenomics/mhcmatch.git@${MHCMATCH_REF}"

# Bootstrap the reference ligand panel from the public HF dataset (both tiers, ~12 MB + ~4 MB) into
# the image's huggingface_hub cache, so from_pmhc() resolves it offline at runtime.
RUN python -c "from mhcmatch.store import fetch_pmhc; fetch_pmhc('full'); fetch_pmhc('shortlist')"

# Sanity: the CLI and panel resolve at build time (fails the build early if the ref/panel is wrong).
RUN mhcmatch --help >/dev/null && \
    python -c "from mhcmatch import Store; Store.from_pmhc(tier='shortlist', species='human', classes=('mhc1',))"

ENTRYPOINT ["mhcmatch"]
CMD ["--help"]
