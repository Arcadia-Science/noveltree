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
        cmake \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        make \
        python3

# Clone the VeryFastTree repo build, and install
RUN git clone -b v3.2.1 https://github.com/citiususc/veryfasttree.git && \
    cd veryfasttree && \
    cmake . && \
    make
    
# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the VeryFastTree binary to path
ENV PATH="/veryfasttree/:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
