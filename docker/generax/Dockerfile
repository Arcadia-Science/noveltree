# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        bison \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        flex \
        git \
        libopenmpi-dev \
        libcurl4-openssl-dev \
        libgmp3-dev \
        libssl-dev \
        openmpi-bin \
        openmpi-doc

# Clone the GeneRax repository and install with MPI
RUN git clone --recursive https://github.com/BenoitMorel/GeneRax && \
    cd GeneRax/ && \
    git checkout 19604b7115c8ff6306f04d69572483de195600a8 && \
    bash ./install.sh

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the GeneRax binary to path
ENV PATH="/GeneRax/build/bin:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
