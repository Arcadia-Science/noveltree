# RAAS Code Vendored from raas-organism-prioritization

This directory contains code vendored from the [raas-organism-prioritization](https://github.com/Arcadia-Science/raas-organism-prioritization) repository for use in the Noveltree pipeline.

## Source Information

- **Repository:** https://github.com/Arcadia-Science/raas-organism-prioritization
- **Commit:** `7619fab7fa47e5f9da2dd2157b76cac6eb125ec0`
- **Date:** 2024-11-14 14:42:42 -0800
- **Commit Message:** Merge pull request #27 from Arcadia-Science/ap/aa-stat-calc

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
- `phylo_correction.cpp` - Phylogenetic GLS transformation (Cholesky decomposition)
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

## Noveltree Pipeline Differences from RAAS

### Output Aggregation

The Noveltree pipeline aggregates per-gene-family outputs into `results/phylo_dist_aggregated/all_protein_comparisons.tsv.gz`.

**Key difference from RAAS:**
- RAAS filters aggregated results to disease-associated genes only (requires ClinVar data)
- Noveltree includes ALL protein comparisons without filtering (no ClinVar integration)
- Disease annotation columns (`disease_name`, etc.) are `NA` since pipeline uses `clinvar = NULL`

**To add disease filtering:** Integrate ClinVar data by updating `modules/local/phylo_dist.nf` to pass `clinvar` parameter instead of NULL, then filter the aggregated output to `!is.na(disease_name)`.
