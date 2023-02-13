# Use an Ubuntu base image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libopenmpi-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        openmpi-bin \
        openmpi-doc

# Clone the GeneRax repository and install with MPI
RUN git clone --recursive https://github.com/BenoitMorel/Asteroid && \
    cd Asteroid/ && \
    git checkout 3aae117df3353c28c6a07d58a4c8b0ab290f974f && \
    bash ./install.sh

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the Asteroid binary to path
ENV PATH="/Asteroid/build/bin:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
