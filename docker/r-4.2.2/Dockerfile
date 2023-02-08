# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent R from asking for your geographic area
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        curl \
        build-essential \
        software-properties-common && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update && \
    apt-get install -y \
        r-base &&

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
