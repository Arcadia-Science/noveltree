require(phytools)
require(Rcpp)
require(RcppArmadillo)
Rcpp::sourceCpp("phylo_correction.cpp")
Rcpp::sourceCpp("compute_inverse_covariance.cpp")
Rcpp::sourceCpp("compute_chunk_distances.cpp")
# The following R-script contains a number of scripts that have been written
# to, as efficiently as possible, calculate multivariate distances between
# species or proteins, given a set of trait measurements and a corresponding
# species tree or gene family tree.

# Numerous functions have been written in Rcpp to take advantage of C++'s
# efficiency, particularly for the mahalanobis distance calculations among
# thousands of proteins and thousands of gene families.

# These functions get to be a bit complex, so I've attemped to walk through
# what they are, and how they work together below. The two core functions are:
# 1) "phylo_gls_transform":
#    - This function carries out the phylogenetic GLS transformation described
#      in Butler et al., 2000 (described in more detail below). In brief,
#      this function takes trait data and a corresponding phylogeny and returns
#      transformed data that represent the residual trait variation after
#      accounting for variance explained by phylogeny.
#    - NOTE - this is effectively a wrapper function around the Rcpp
#             function "phylo_correction"
# 2) "pairwise_mahalanobis":
#    - This takes the transformed data and the tree, and calculates the
#      pairwise mahalanobis distances. Distances are calculated using data
#      corresponding to terminal nodes (e.g. species) by default. Otherwise
#      distances can be calculated using data corresponding to all terminal
#      AND internal nodes by specifying include_internal = FALSE. Doing so
#      is appropriate when calculating distances among gene families using
#      gene family evolutionary event counts.
#    - NOTE - Internally this calls a number of internal functions (described
#             in more detail later):
#             - "compute_inverse_covariance"
#             - "get_dist_pair_chunks"
#             - "process_all_chunks"
# 3) "phylo_correction":
#    - This Rcpp function takes as input the phylogenetic covariance matrix and
#      traits, returning either a vector or matrix of phylogenetic GLS
#      transformed data if one or multiple traits are provided respectively
# 4) "compute_inverse_covariance":
#    - This Rcpp function calculates the inverse of the covariance matrix of
#      transformed data for use in calculation of the mahalanobis distances.
#      It is implemented in Rcpp for the increased efficiency of matrix algebra
# 5) "get_dist_pair_chunks":
#    - Helper function to facilitate the calculation of mahalanobis distances.
#    - Generates all pairwise comparisons (among indices), and then splits
#      these pairs into smaller groups (chunks) of a specified size.
#      Each chunk is a matrix of pairs, and the function returns a list of
#      these chunk matrices.
# 6) "process_all_chunks":
#    - Takes the list of "chunked" pairwise comparisons produced by
#      "get_dist_pair_chunks" and the inverse phenotypic covariance matrix
#      produced by "compute_inverse_covariance", and calculates the pairwise
#      mahalanobis distances in parallel using mclapply returning the full
#      matrix.
#    - Internally, calls the Rcpp function "compute_chunk_distances" that
#      calculates the actual distances for each chunk.
# 7) "compute_chunk_distances":
#     - This function calculates the actual pairwise mahalanobis distances
#       using the matrix of data, inverse of the phenotypic covariance matrix,
#       and matrix of indices indicating the pair of observations we want to
#       calculate distances among (i.e. the chunk)

# Define a function that conducts the phylogenetic GLS transformation,
# returning the residual trait variation after accounting for the phylogenetic
# effect for a set of samples that correspond only to the tips of the tree (if
# include_internal == FALSE), or the tips of the tree and internal nodes
# (i.e. if conducting phylogenetic profiling using counts of gene
# duplication / transfer, or loss inferred using GeneRax)
phylo_gls_transform <- function(data_matrix, tree, include_internal = FALSE) {
  # data_matrix: Assumed to be of class "matrix", containing
  # the phenotypes/traits you would like to transform, with dimnames
  # corresponding to trait and node names - the function checks whether node
  # names are in columns or rows and handles internally.
  # tree: Assumed to be a tree of class "phylo" with nodes
  # labeled corresponding to the tips and internal nodes
  # if include_internal == TRUE
  # include_internal: determines whether the traits include
  # internal nodes, or if you would like to return the transformed
  # data for all (non-root nodes)

  # Obtain the phylogenetic VCV for species and (if selected) internal nodes:
  phylo_vcv <- phytools::vcvPhylo(tree, anc.nodes = include_internal)

  # Identify if node ids are in columns this will be used to return
  # traits in the same format as they were input to the function
  nodes_as_columns <- all(tree$tip.label %in% colnames(data_matrix))

  if (include_internal == TRUE) {
    colnames(phylo_vcv) <- c(tree$tip.label, tree$node.label[-1])
    rownames(phylo_vcv) <- c(tree$tip.label, tree$node.label[-1])

    # Identify the root node:
    root_node_name <- tree$node.label[1]

    # And prepare the data for transformation:
    data_matrix <- prep_data(data_matrix, tree,
                             include_internal, root_node_name)
  } else {
    data_matrix <- prep_data(data_matrix, tree)
  }

  # Now, obtain the phylogenetic GLS transformed data
  transf_data <-
    do.call(cbind, phylo_correction(data_matrix, phylo_vcv)) # nolint
  dimnames(transf_data) <- list(rownames(data_matrix), colnames(data_matrix))

  # Return the phylogenetically transformed data in the original format
  if (nodes_as_columns) {
    transf_data <- t(transf_data)
  }
  return(transf_data)
}


# A separate function to do the same as above, but given a list of
# input data matrices
phylo_gls_transform_list <-
  function(data_matrix_list, tree, include_internal = FALSE) {
    # data_matrix_list: A list of traits, each stored as class "matrix",
    # containing traits corresponding to either terminal or both terminal
    # and internal nodes of the species tree. All input matrices must be
    # structured in the same way.
    # tree: Assumed to be a species tree of class "phylo"
    # with nodes labeled corresponding to the species within the
    # count_matrices.
    # include_internal: determines whether the count data include
    # internal nodes, or if you would like to return the transformed
    # data for all (non-root nodes)

    # Obtain the phylogenetic VCV for species and (if selected) internal nodes:
    phylo_vcv <- phytools::vcvPhylo(tree, anc.nodes = include_internal)

    # Identify if node ids are in columns this will be used to return
    # traits in the same format as they were input to the function
    # Again, assume all matrices are structured in the same way
    nodes_as_columns <-
      all(tree$tip.label %in% colnames(data_matrix_list[[1]]))

    if (include_internal == TRUE) {
      colnames(phylo_vcv) <- c(tree$tip.label, tree$node.label[-1])
      rownames(phylo_vcv) <- c(tree$tip.label, tree$node.label[-1])

      # Identify the root node:
      root_node_name <- tree$node.label[1]

      # And prepare the data for transformation:
      data_matrix_list <-
        lapply(data_matrix_list, function(x) {
          prep_data(x, tree, include_internal, root_node_name)
        })
    } else {
      data_matrix_list <-
        lapply(data_matrix_list,
               function(x) prep_data(x, tree))
    }

    # Now, obtain the phylogenetic GLS transformed data for each
    # matrix in the list
    transf_data <- list()
    for (i in seq_along(data_matrix_list)) {
      transf_data[[i]] <-
        do.call(cbind, phylo_correction(data_matrix_list[[i]], # nolint
                                        phylo_vcv))
      dimnames(transf_data[[i]]) <-
        list(rownames(data_matrix_list[[i]]), colnames(data_matrix_list[[i]]))
    }

    # Now combine each transformed matrix into a single matrix, again,
    # maintaining the original format
    if (nodes_as_columns) {
      # Return the phylogenetically transformed data in the original format
      transf_data <- lapply(transf_data, t)
      transf_data <- do.call(cbind, transf_data)
    } else {
      transf_data <- do.call(rbind, transf_data)
    }
    return(transf_data)
  }

# A function to prepare trait data for the phylo_correction function
prep_data <- function(data_matrix, tree,
                      include_internal = FALSE,
                      root_node_name = NULL) {
  # data_matrix: single matrix of traits corresponding to species tree
  # tree: phylogeny as class phylo, with nodes labeled to match colnames
  # in data_matrix
  # Ensure our data is a matrix:
  if (!is.matrix(data_matrix)) {
    data_matrix <- as.matrix(data_matrix)
  }
  # Determine whether tips are stored in columns or rows
  if (all(tree$tip.label %in% colnames(data_matrix))) {
    # transpose data matrix so species are stored in rows
    data_matrix <- t(data_matrix)
  }

  if (include_internal == TRUE) {
    # Remove the root node from the data (if included):
    data_matrix <-
      data_matrix[-which(rownames(data_matrix) == root_node_name), ]
    # Ensure the data matches the order of the species tree
    data_matrix <-
      data_matrix[match(c(tree$tip.label, tree$node.label[-1]),
                        rownames(data_matrix)), ]
  } else {
    # Ensure the data matches the order of the species tree
    data_matrix <- data_matrix[match(tree$tip.label, rownames(data_matrix)), ]
  }
  return(data_matrix)
}

# A function to calculate all pairwise Mahalanobis distances
pairwise_mahalanobis <- function(data_matrix) {
  # data_matrix: for now, output data (class matrix) from phylo_gls_transform
  # function, but on principle can be applied to any any properly formatted
  # data
  # NOTE: data_matrix must be formatted such that observations among which
  # distances are calculated must be stored by-row.
  # Calculate the inverse of the covariance matrix for calculation of #
  # mahalanobis distances
  inv_cov <- compute_inverse_covariance(data_matrix) # nolint

  # Break generate a list of all possible unique pairwise comparisons,
  # dividing them into chunks of 10000
  pair_chunks <- get_dist_pair_chunks(nrow(data_matrix), 10000)

  # Now, calculate distances for each chunk with Rcpp, parallelizing across
  # chunks and returning a symmetrical distance matrix
  distance_matrix <- process_all_chunks(pair_chunks, data_matrix, inv_cov)
  rownames(distance_matrix) <- rownames(data_matrix)
  colnames(distance_matrix) <- rownames(data_matrix)
  return(distance_matrix)
}

# A function to facilitate the parallelization of mahalanobis distances
# by breaking the pairwise sets of comparisons into a list of "chunks"
get_dist_pair_chunks <- function(n, chunk_size) {
  total_pairs <- combn(1:n, 2)
  num_chunks <- ceiling(ncol(total_pairs) / chunk_size)
  split_cols <-
    split(seq_len(ncol(total_pairs)),
          rep(1:num_chunks, each = chunk_size, len = ncol(total_pairs)))
  lapply(split_cols, function(cols) total_pairs[, cols, drop = FALSE])
}

# Calculate multivariate distances among samples using the indices from a
# single "chunk"
process_single_chunk <- function(pair_chunk, data_matrix, inv_cov) {
  chunk_distances <-
    compute_chunk_distances(data_matrix, inv_cov, pair_chunk) # nolint
  list(pair_chunk = pair_chunk, distances = chunk_distances)
}

# A function that takes the list of "chunked" pairwise comparisons produced by
# `get_dist_pair_chunks` and the inverse phenotypic covariance matrix produced
# by `compute_inverse_covariance`, and calculates the pairwise mahalanobis
# distances in parallel using mclapply returning the full matrix.
process_all_chunks <- function(pair_chunks, data_matrix, inv_cov) {
  n <- nrow(data_matrix)
  distance_matrix <- matrix(NA, n, n)

  # Process the chunks in parallel for efficiency
  chunk_results <-
    parallel::mclapply(pair_chunks, process_single_chunk,
                       data_matrix = data_matrix, inv_cov = inv_cov,
                       mc.cores = parallel::detectCores())

  # Populate a distance matrix and set the diagonal equal to zero
  for (result in chunk_results) {
    pair_chunk <- result$pair_chunk
    distances <- result$distances
    row_indices <- pair_chunk[1, ]
    col_indices <- pair_chunk[2, ]
    for (i in seq_along(row_indices)) {
      distance_matrix[row_indices[i], col_indices[i]] <- distances[i]
      distance_matrix[col_indices[i], row_indices[i]] <- distances[i]
    }
  }
  diag(distance_matrix) <- 0
  return(distance_matrix)
}
