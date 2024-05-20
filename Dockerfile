# Use Ubuntu 20.04 as the base image for better ARM support
FROM ubuntu:20.04

LABEL maintainer="Your Name <your.email@example.com>"

ENV DEBIAN_FRONTEND=noninteractive

# Install required system packages
RUN apt-get update -y \
    && apt-get install -y wget build-essential zlib1g-dev git

# Download and install Miniforge (a Miniconda variant for ARM)
RUN wget -O /tmp/miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh \
    && chmod +x /tmp/miniforge.sh \
    && /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh

# Add Conda to the PATH
ENV PATH /opt/conda/bin:$PATH

# Copy the environment file into the image
COPY environment.yml /environment.yml

# Create Conda environment and install dependencies
RUN conda env create --name bioinfo_env --file /environment.yml \
    && conda clean -a

# Add conda environment to the PATH
ENV PATH /opt/conda/envs/bioinfo_env/bin:$PATH
ENV CONDA_DEFAULT_ENV=bioinfo_env
ENV CONDA_PREFIX=/opt/conda/envs/bioinfo_env

# Copy the R package to the Docker image
COPY . /usr/local/lib/R/site-library/barbac

# Set the working directory
WORKDIR /data

# Set the entrypoint to Rscript
ENTRYPOINT ["Rscript"]

# Define the default command
CMD ["-e", "barbac::qc_merge"]