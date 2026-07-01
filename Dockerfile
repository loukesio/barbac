# -----------------------------------------------------------------------------
# barbac Docker image
#
# Reproducible environment for DNA barcode lineage analysis.
# Contains:
#   * R + Bioconductor (via bioconductor/bioconductor_docker)
#   * CLI pipeline tools: FastQC, PEAR, minimap2, samtools, MultiQC
#   * The `barbac` R package installed from the working copy in this repo
#
# Usage (once published to GHCR):
#   docker pull ghcr.io/loukesio/barbac:latest
#   docker run --rm -it -v $(pwd):/data ghcr.io/loukesio/barbac:latest R
#
# Build locally:
#   docker build -t barbac:dev .
# -----------------------------------------------------------------------------
FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="barbac" \
      org.opencontainers.image.description="End-to-end DNA barcode lineage analysis pipeline (R + FASTQ tools)" \
      org.opencontainers.image.source="https://github.com/loukesio/barbac" \
      org.opencontainers.image.licenses="GPL-2.0-or-later"

# -----------------------------------------------------------------------------
# System-level CLI tools available in Debian/Ubuntu repos
# -----------------------------------------------------------------------------
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        fastqc \
        samtools \
        minimap2 \
        curl \
        bzip2 \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# -----------------------------------------------------------------------------
# PEAR and MultiQC via bioconda through micromamba. PEAR isn't in Debian
# repos (direct download requires registration); the newer Debian base
# used by bioconductor_docker blocks system-Python `pip install` under
# PEP 668, so multiqc goes through conda too.
# -----------------------------------------------------------------------------
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/envs/tools/bin:${MAMBA_ROOT_PREFIX}/bin:${PATH}
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C /usr/local/bin bin/micromamba && \
    micromamba create -y -n tools -c bioconda -c conda-forge \
        pear=0.9.6 \
        multiqc=1.25.1 && \
    micromamba clean -a -y

# -----------------------------------------------------------------------------
# barbac R package (deps + source install + smoke test)
# -----------------------------------------------------------------------------
WORKDIR /pkg
COPY . /pkg

RUN R --no-save <<'RSCRIPT'
options(
  repos       = BiocManager::repositories(),
  Ncpus       = parallel::detectCores(),
  install.packages.check.source = "no"
)
install.packages(c("remotes", "devtools"))
remotes::install_deps(dependencies = TRUE, upgrade = "never")
# ltc is a github-only palette package used by the vignette examples
remotes::install_github("loukesio/ltc_palettes", upgrade = "never")
devtools::install(quiet = TRUE, upgrade = "never")

# Smoke-test: fail the build if barbac cannot load or export its API
stopifnot(
  is.function(barbac::super_cluster2),
  is.function(barbac::barbac_ts_area),
  is.function(barbac::run_cli_pipeline)
)
message("barbac loaded successfully.")
RSCRIPT

# Sensible working directory for user data mounts
WORKDIR /data

CMD ["R"]
