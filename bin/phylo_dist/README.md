# Organism Selection Code Vendored from 2024-organismal-selection

This directory contains code vendored from the [2024-organismal-selection](https://github.com/Arcadia-Science/2024-organismal-selection) repository for use in the Noveltree pipeline.

## Source Information

- **Repository:** https://github.com/Arcadia-Science/2024-organismal-selection

## Purpose

These scripts provide two key functionalities for the Noveltree pipeline:

1. **Physicochemical Property Calculation** (`genefam_aa_summaries.py`)

   - Calculates amino acid composition and physicochemical properties for protein families
   - Used by the `PHYSICOCHEMICAL_PROPS` module

2. **Phylogenetically-Corrected Protein Distance Analysis** (R and C++ files)
   - Performs phylogenetic GLS transformation to remove phylogenetic signal
   - Calculates pairwise Mahalanobis distances between proteins
   - Conducts statistical tests for protein similarity
   - Used by the `PHYLO_DIST` module

## Files

### Python Scripts

- `genefam_aa_summaries.py` - Calculate physicochemical properties from MSA files

### R Scripts

- `phylo_multivariate_distance_functions.R` - Core phylogenetic GLS transformation and distance calculation
- `protein_distance_calculation_functions.R` - Main functions for protein distance analysis
- `protein_dist_permutation_tests.R` - Permutation testing for statistical significance
- `simple_protein_dist_signif_tests.R` - Z-score and Wilcoxon tests for protein similarity

### C++ Files (used via Rcpp)

- `phylo_correction.cpp` - Phylogenetic GLS transformation
- `calculate_dist_stats.cpp` - Calculate summary statistics for distance matrices
- `compute_inverse_covariance.cpp` - Compute inverse covariance for Mahalanobis distances
- `compute_chunk_distances.cpp` - Pairwise Mahalanobis distance calculation
- `perm_test_across_non_ref.cpp` - Permutation test across non-reference proteins
- `perm_test_within_non_ref.cpp` - Permutation test within non-reference proteins

## Modifications

- **None** - These are exact copies from the source repository at the specified commit.

## Dependencies

### R Packages

- Rcpp, RcppArmadillo
- dplyr, tidyr
- ape, phytools, geiger

### Python Packages

- biopython
- pandas
- argparse (standard library)

## Integration Notes

- All R scripts use relative paths to source C++ files via `Rcpp::sourceCpp()`
- C++ files are compiled at Docker build time
- Scripts are designed to work from this directory location

---
