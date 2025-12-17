# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND=noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        python3.9 \
        python3.9-distutils \
        wget && \
    wget https://bootstrap.pypa.io/get-pip.py && \
    python3.9 get-pip.py && \
    rm get-pip.py && \
    python3.9 -m pip install clipkit==2.1.1 seqmagick==0.8.4 && \
    ln -sf /usr/bin/python3.9 /usr/bin/python3 && \
    ln -sf /usr/local/bin/pip3.9 /usr/bin/pip3

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
