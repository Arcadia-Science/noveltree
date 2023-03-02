# Use an Ubuntu base image
FROM ubuntu:20.04

# Set path for miniconda
ENV PATH=/opt/miniconda3/bin:$PATH

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        python3 \
        wget

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

# Install veryfasttree
RUN conda install --quiet -c bioconda veryfasttree

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
