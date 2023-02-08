# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent R from asking for your geographic area
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        bison \
        flex \
        ca-certificates \
        cmake \
        openmpi-bin \
        openmpi-doc \
        libopenmpi-dev \
        build-essential \
        libcurl4-openssl-dev \
        libgmp3-dev \
        libssl-dev \
        git
        
RUN git clone --recursive https://github.com/BenoitMorel/GeneRax && \
    cd GeneRax/ && \
    bash ./install.sh

ENV PATH="/GeneRax/build/bin:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
