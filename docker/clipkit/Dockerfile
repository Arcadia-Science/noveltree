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
        python3.9 \
        python3-pip \
        wget && \
    pip3 install clipkit==2.1.1 biopython==1.79 seqmagick==0.8.4

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
