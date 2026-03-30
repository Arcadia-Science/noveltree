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

# Install FAMSA from GitHub
RUN git clone https://github.com/refresh-bio/FAMSA.git /opt/FAMSA && \
    cd /opt/FAMSA && \
    make && \
    mv bin/famsa /usr/local/bin/ && \
    rm -rf /opt/FAMSA

# Set up working directory
WORKDIR /data

# Default command
CMD ["/bin/bash", "-l"]
