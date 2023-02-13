FROM arcadiascience/rbase:4.2.2

# Install the required dependencies and R
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        dirmngr \
        gpg-agent \
        libcurl4-openssl-dev \
        libicu-dev \
        libssl-dev \
        libtiff-dev \
        libxml2-dev \
        make \
        pandoc \
        zlib1g-dev

# Install R packages required to run cogeqc module
RUN Rscript -e "install.packages('tidyverse', version = '1.3.2', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('BiocManager', version = '3.16', repos='https://cloud.r-project.org/')"
RUN Rscript -e "BiocManager::install('cogeqc')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
