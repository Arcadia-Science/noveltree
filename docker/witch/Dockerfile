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

# Build and install
RUN cd WITCH && \ 
    pip install -r requirements.txt && \
    python setup.py
    
# Modify the backbone construction
COPY ./modified_witch_backbone.py .
# Retain the original script, with a slightly modified name
RUN cp /WITCH/gcmm/backbone.py /WITCH/gcmm/backbone_original.py
RUN mv modified_witch_backbone.py /WITCH/gcmm/backbone.py

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
