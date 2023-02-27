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
        git

RUN git clone https://github.com/KatyBrown/CIAlign.git && \
    cd CIAlign/ && \
    sed "s/pandas==/pandas>=/g" requirements.txt

RUN pip3 install Cython && \
    python3.9 pip install -r requirements.txt

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the CIAlign binary directory to path
ENV PATH="CIAlign/:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
