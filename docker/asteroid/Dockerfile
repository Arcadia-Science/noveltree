FROM arcadiascience/rbase_4.2.2:0.0.1

# Prevent interactive querying of geographic area during build
ENV DEBIAN_FRONTEND noninteractive

# Install the required dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential \
        ca-certificates \
        cmake \
        git \
        libopenmpi-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        openmpi-bin \
        openmpi-doc \
        python3.7 \
        python-is-python3 \
        pip

# Clone the Asteroid repository and install
RUN git clone --recursive https://github.com/BenoitMorel/Asteroid && \
    cd Asteroid/ && \
    git checkout 3aae117df3353c28c6a07d58a4c8b0ab290f974f && \
    bash ./install.sh
    
# Clone the DISCO repository and install, along with dependencies
RUN git clone https://github.com/JSdoubleL/DISCO.git && \
    cd DISCO/ && \
    git checkout 20e10c33296a870a3f60c695afced69dcf36489b && \
    pip install treeswift

# Install ape and phytools for species tree rooting
RUN Rscript -e "install.packages('ape', version = '5.7.1', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('phytools', version = '1.5.1', repos='https://cloud.r-project.org/')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the Asteroid binary to path
ENV PATH="/Asteroid/build/bin:$PATH"

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
