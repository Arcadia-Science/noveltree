# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        python3.8 \
        python3-pip \
        python-is-python3 \
        git

RUN git clone -b v1.1.0 https://github.com/KatyBrown/CIAlign.git && \
    cd CIAlign/ && \
    pip3 install -r requirements.txt && \
    python setup.py build && \
    python setup.py install 
    
# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the CIAlign binarycd directory to path
ENV PATH="CIAlign/build/scripts-3.8:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
