# Use an Ubuntu base image
FROM ubuntu:20.04

# Prevent R from asking for your geographic area
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies and R
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        build-essential \
        curl \
        dirmngr \
        gpg-agent \
        libcurl4-openssl-dev \
        libicu-dev \
        libssl-dev \
        libtiff-dev \
        libxml2-dev \
        make \
        pandoc \
        software-properties-common \
        zlib1g-dev 

# Install R 4.2.2
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get update && \
    apt-get install -y \
        r-base

# Install R packages required to run cogeqc module
RUN Rscript -e "install.packages(c('tidyverse', 'BiocManager'), repos='https://cloud.r-project.org/')"
RUN Rscript -e "BiocManager::install('cogeqc')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
