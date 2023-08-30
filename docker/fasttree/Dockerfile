FROM arcadiascience/rbase_4.2.2:1.0.0

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        gcc \
        git \
        libopenmpi-dev \
        libcurl4-openssl-dev \
        libgmp3-dev \
        libssl-dev \
        openmpi-bin \
        openmpi-doc \
        wget

# Download the FastTree source code and compile with OpenMP for
# parallel execution and use double precision to prevent rounding/
# numerical issues
RUN wget http://microbesonline.org/fasttree/FastTree.c && \
    gcc -DUSE_DOUBLE -DOPENMP -fopenmp -O3 -finline-functions -funroll-loops -Wall -o FastTreeDblMP FastTree.c -lm && \
    chmod u+x FastTreeDblMP && \
    cp FastTreeDblMP /bin/

# Install ape to resolve polytomies
RUN Rscript -e "install.packages('ape', version = '5.6.2', repos='https://cloud.r-project.org/')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
