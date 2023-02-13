FROM arcadiascience/rbase:4.2.2

# Install the required dependencies and R
RUN apt-get update && \
    apt-get install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        zlib1g-dev

RUN Rscript -e "install.packages('BiocManager', version = '3.16', repos='https://cloud.r-project.org/')"
RUN Rscript -e "BiocManager::install('UniProt.ws')"

RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
