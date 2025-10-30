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
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages(c('Rcpp', 'RcppArmadillo', 'dplyr', 'tidyr', 'ape', 'phytools', 'geiger'), repos='https://cloud.r-project.org/')"

# Create directory for source scripts
RUN mkdir -p /opt/phylo_dist

# Copy src R and C++ scripts from the repository root
# Note: Docker build context is the repository root, so we reference bin/phylo_dist/ from there
COPY ../../bin/phylo_dist/*.R /opt/phylo_dist/
COPY ../../bin/phylo_dist/*.cpp /opt/phylo_dist/

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
