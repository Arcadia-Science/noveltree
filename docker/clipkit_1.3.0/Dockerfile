# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent R from asking for your geographic area
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        python3 \
        python3-pip \
        wget && \
    pip3 install clipkit

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
