FROM arcadiascience/rbase_4.2.2:1.0.0

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
    zlib1g-dev \
    python3-pip

RUN pip install biopython

# Install R package required to run protein_annotation module
RUN Rscript -e "install.packages('queryup', version = '1.0.5', repos='https://cloud.r-project.org/')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
