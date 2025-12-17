FROM rocker/r-ver:4.3.1

# Install system dependencies
RUN apt-get update && \
    apt-get install -y \
    procps \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgfortran5 \
    gfortran \
    g++ \
    cmake \
    git \
    wget \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda for conda package management
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -afy

# Add conda to PATH and set CONDA_PREFIX for R scripts
ENV PATH="/opt/conda/bin:${PATH}"
ENV CONDA_PREFIX="/opt/conda"

# Accept conda Terms of Service and configure channels
RUN conda config --set channel_priority flexible && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install treePL via conda
# Note: conda package installs binary to /opt/conda/treePL, need to symlink to bin/
RUN conda install -c genomedk treepl -y && \
    conda clean -afy && \
    ln -s /opt/conda/treePL /opt/conda/bin/treePL

# Download and compile PATHd8 from source
# PATHd8 is a C program for phylogenetic dating, official site: www2.math.su.se/PATHd8/
RUN wget https://www2.math.su.se/PATHd8/PATHd8.zip -O /tmp/PATHd8.zip && \
    unzip /tmp/PATHd8.zip -d /opt/PATHd8_src && \
    cd /opt/PATHd8_src && \
    gcc PATHd8.c -O3 -lm -o PATHd8 && \
    mv PATHd8 /usr/local/bin/ && \
    cd / && \
    rm -rf /opt/PATHd8_src /tmp/PATHd8.zip

# Install R packages
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'dplyr', 'tidyr', 'ape', 'phytools', 'geiger'), repos='https://cloud.r-project.org/')"

# Create directory for source scripts
RUN mkdir -p /opt/phylo_dist

# Copy src R and C++ scripts from the repository root
# Note: Docker build context is the repository root, so we reference bin/phylo_dist/ from there
# This includes time_calibrate_species_tree.R which is now in bin/phylo_dist/
COPY bin/phylo_dist/*.R /opt/phylo_dist/
COPY bin/phylo_dist/*.cpp /opt/phylo_dist/

# Set working directory and compile C++ code
WORKDIR /opt/phylo_dist
RUN R -e "Rcpp::sourceCpp('phylo_correction.cpp')" && \
    R -e "Rcpp::sourceCpp('calculate_dist_stats.cpp')" && \
    R -e "Rcpp::sourceCpp('compute_inverse_covariance.cpp')" && \
    R -e "Rcpp::sourceCpp('compute_chunk_distances.cpp')" && \
    R -e "Rcpp::sourceCpp('perm_test_across_non_ref.cpp')" && \
    R -e "Rcpp::sourceCpp('perm_test_within_non_ref.cpp')"

# Add scripts to PATH
ENV PATH="/opt/phylo_dist:${PATH}"

# Set the default command to launch an interactive shell
CMD ["/bin/bash", "-l"]
