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
        python3.9 \
        python-is-python3 \
        git

# Clone the WITCH github repo
RUN git clone -b v0.3.0 https://github.com/c5shen/WITCH.git

# Install the required Python dependencies
# Here we're not using the requirements.txt at the root of the repo
# because the repo doesn't pin versions, which leads to dependency conflicts
RUN pip install --no-cache-dir \
    configparser==5.3.0 \
    DendroPy==4.5.2 \
    numpy==1.24.2 \
    psutil==5.9.4

# Build WITCH
RUN cd WITCH && \
    python setup.py

# Retain the original backbone construction script:
RUN mv /WITCH/gcmm/backbone.py /WITCH/gcmm/backbone_original.py

# Modify the backbone construction
COPY modified_witch_backbone.py /WITCH/gcmm/backbone.py

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
