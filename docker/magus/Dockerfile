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
        git \
        pip \
        python3.8 \
        python-is-python3 \
        wget
        
# Clone MAGUS and VeryFastTree github repos
RUN git clone https://github.com/Arcadia-Science/MAGUS-VFT.git

# Build and install
RUN cd MAGUS-VFT && \
    git checkout 435d142cb196c9ffca01a937e9751ce3ab3417f7 && \
    python setup.py build && \
    python setup.py install
    
# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
