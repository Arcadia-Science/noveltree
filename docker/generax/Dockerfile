FROM arcadiascience/rbase_4.2.2:0.0.1

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        bison \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        flex \
        git \
        libopenmpi-dev \
        libcurl4-openssl-dev \
        libgmp3-dev \
        libssl-dev \
        openmpi-bin \
        openmpi-doc

# Clone the GeneRax repository and install with MPI
RUN git clone --recursive https://github.com/BenoitMorel/GeneRax && \
    cd GeneRax/ && \
    git checkout 5d52c7d08f8ad15f49b1d85507773602e46e93a8 && \
    bash ./install.sh

# Install ape and phytools for species tree rooting
RUN Rscript -e "install.packages('ape', version = '5.7.1', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('phytools', version = '1.5.1', repos='https://cloud.r-project.org/')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the GeneRax binary to path
ENV PATH="/GeneRax/build/bin:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
