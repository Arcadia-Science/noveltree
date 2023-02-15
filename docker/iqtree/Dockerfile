# Use an Ubuntu base image
FROM ubuntu:20.04

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential \
        ca-certificates \
        cmake \
        wget \
        libcurl4-openssl-dev \
        libeigen3-dev \
        libssl-dev

# Clone the GeneRax repository and install with MPI
RUN wget https://github.com/iqtree/iqtree2/archive/refs/tags/v2.2.0.5.tar.gz \
    tar -xvf v2.2.0.5.tar.gz && \
    rm v2.2.0.5.tar.gz && \
    cd iqtree2-2.2.0.5/ && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
