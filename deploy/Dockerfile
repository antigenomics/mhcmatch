# mhcmatch epitope-prediction image for the Gamaleya nextflow_vaccine pipeline.
#
# Packages `mhcmatch predict` + its seqtree C++ core + the reference ligand panel, so a process can
# run:  mhcmatch predict <fasta> --alleles <...> --cls <mhc1|mhc2> --scored-csv <out> --native <out>
# with no host data mounts (the panel is baked; $MHCMATCH_PMHC points at it).
#
# Build (from a context that also holds the two panel files — see deploy/README.md):
#   cp ~/hf/pmhc_data/pmhc/pmhc_full.tsv.gz ~/hf/pmhc_data/pmhc/pmhc_shortlist.tsv.gz deploy/
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

# Bake the reference ligand panel (isalgo/pmhc_data). Small enough to ship (~12 MB + ~4 MB).
# from_pmhc(None, tier=T) resolves $MHCMATCH_PMHC/pmhc_<T>.tsv.gz, so drop the files here directly.
COPY pmhc_full.tsv.gz pmhc_shortlist.tsv.gz /opt/pmhc_data/
ENV MHCMATCH_PMHC=/opt/pmhc_data

# Sanity: the CLI and panel resolve at build time (fails the build early if the ref/panel is wrong).
RUN mhcmatch --help >/dev/null && \
    python -c "from mhcmatch import Store; Store.from_pmhc(tier='shortlist', species='human', classes=('mhc1',))"

ENTRYPOINT ["mhcmatch"]
CMD ["--help"]
