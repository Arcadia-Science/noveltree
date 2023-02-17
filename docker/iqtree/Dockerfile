# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential \
        ca-certificates \
        clang-tidy \
        cmake \
        git \
        wget \
        libboost-all-dev \
        libcurl4-openssl-dev \
        libeigen3-dev \
        libssl-dev
        

# Clone the GeneRax repository and install with MPI
RUN git clone https://github.com/iqtree/iqtree2.git --branch v2.2.0.5 && \
    cd iqtree2 && \
    mkdir build && \
    cd build && \
    cmake -DUSE_LSD2=OFF .. && \
    make && \
    make install

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
