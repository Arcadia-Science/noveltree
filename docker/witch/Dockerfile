# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND=noninteractive

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        cmake \
        libcurl4-openssl-dev \
        libssl-dev \
        git \
        pip \
        python3.9 \
        python-is-python3

# Install WITCH v1.0.10 via pip (includes MAGUS backbone aligner and FastTree2)
RUN pip install --no-cache-dir witch-msa==1.0.10

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
