# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        cmake \
        libcurl4-openssl-dev \
        libssl-dev \
        python3 \
        python3-pip \
        wget \
        git

# Clone MAGUS and VeryFastTree github repos
RUN git clone -b ap/add-vft https://github.com/Arcadia-Science/MAGUS-VFT

# Build and install
RUN cd MAGUS-VFT && \
    python3 setup.py build && \
    python3 setup.py install
    
# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
