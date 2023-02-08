# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent R from asking for your geographic area
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        make \
        curl \
        pandoc \
        dirmngr \
        gpg-agent

RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        libcurl4-openssl-dev \
        libssl-dev \
        libicu-dev \
        libxml2-dev \
        libtiff-dev \
        zlib1g-dev \
        build-essential \
        software-properties-common

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update && \
    apt-get install -y \
        r-base

RUN Rscript -e "install.packages(c('tidyverse', 'BiocManager'), repos='https://cloud.r-project.org/')"
RUN Rscript -e "BiocManager::install('cogeqc')"

RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
