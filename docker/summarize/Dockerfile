FROM arcadiascience/rbase_4.2.2:1.0.0

ENV MAKEFLAGS="-j4"

# Fixes an issue with timedatectl on Docker.
ENV TZ UTC
RUN echo "TZ=$( cat /etc/timezone )" >> /etc/R/Renviron.site

# System dependencies for the various R packages.
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        build-essential \
        ca-certificates \
        clang-tidy \
        cmake \
        dirmngr \
        gpg-agent \
        libboost-all-dev \
        libcurl4-openssl-dev \
        libeigen3-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libgmp3-dev \
        libharfbuzz-dev \
        libicu-dev \
        libopenmpi-dev \
        libssl-dev \
        libtiff-dev \
        libxml2-dev \
        make \
        pandoc \
        zlib1g-dev

# R packages.
RUN Rscript -e "install.packages('tidyverse', version = '1.3.2')"
RUN Rscript -e "install.packages('BiocManager', version = '3.16')"
RUN Rscript -e "BiocManager::install('cogeqc')"
RUN Rscript -e "install.packages('ape', version = '5.7.1')"
RUN Rscript -e "install.packages('phytools', version = '1.5.1')"
RUN Rscript -e "install.packages('ggtree', version = '3.6.2')"
RUN Rscript -e "install.packages('ggnewscale', version = '0.4.8')"
RUN Rscript -e "install.packages('parallel', version = '4.2.1')"
RUN Rscript -e "install.packages('ggforce', version = '0.4.1')"
RUN Rscript -e "install.packages('plotly', version = '4.10.1')"
RUN Rscript -e "install.packages('ape', version = '5.7.1')"
RUN Rscript -e "install.packages('cowplot', version = '1.1.1')"
RUN Rscript -e "install.packages('reshape2', version = '1.4.4')"
RUN Rscript -e "install.packages('phytools', version = '1.5.1')"
RUN Rscript -e "install.packages('optparse', version = '1.7.4')"
RUN Rscript -e "install.packages('rmarkdown', version = '2.25')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
