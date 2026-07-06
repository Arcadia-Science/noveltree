# Use an Ubuntu base image with newer GCC
FROM ubuntu:22.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
RUN apt-get update && \
    apt-get install -y \
    g++ \
    make \
    cmake \
    git \
    procps \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install FAMSA from GitHub, pinned to the exact commit that ships in the
# published image (FAMSA 2.5.0). The image is named `famsa_2.0.0` for historical
# reasons; the actual tool version built here is 2.5.0 (there is no v2.5.0 tag
# upstream, so we pin the commit SHA directly for reproducibility).
RUN git clone https://github.com/refresh-bio/FAMSA.git /opt/FAMSA && \
    cd /opt/FAMSA && \
    git checkout 5a326d5d262aa5533af07eeeaf41822e67604520 && \
    make && \
    mv bin/famsa /usr/local/bin/ && \
    rm -rf /opt/FAMSA

# Set up working directory
WORKDIR /data

# Default command
CMD ["/bin/bash", "-l"]
