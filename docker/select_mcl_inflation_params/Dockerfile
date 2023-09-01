FROM arcadiascience/rbase_4.2.2:1.0.0

# Install the required dependencies
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        make \
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
        libfreetype6-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libfontconfig1-dev \
        zlib1g-dev

RUN Rscript -e "install.packages('tidyverse', version = '1.3.2', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('cowplot', version = '1.1.1', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('reshape', version = '0.8.9', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('devtools', version = '2.4.5', repos='https://cloud.r-project.org/')"
RUN Rscript -e "devtools::install_github('ahasverus/elbow@cf58784f3a5c17365a1743453a3212cc07a9f216')"

RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
