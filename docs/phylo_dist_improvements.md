# Proposed improvements to phylogenetic distance analysis

## Background

The `PHYLO_DIST` module computes phylogenetically-corrected protein distances within each gene family. The current pipeline:

1. Calculates physicochemical properties for each protein (molecular weight, aromaticity, instability, flexibility, GRAVY, isoelectric point, charge at pH 7, helix/sheet fractions, molar extinction coefficient).
2. Time-calibrates gene family trees by congruifying them with the time-calibrated species tree.
3. Applies a phylogenetic GLS transformation (Butler et al., 2000) to remove trait variance explained by shared evolutionary history, using a Brownian motion covariance model.
4. Computes pairwise Mahalanobis distances on the transformed data to obtain phylogenetically-corrected protein distances.
5. Tests whether individual proteins or species are significantly more similar to a reference species than expected by chance.

Two targeted changes would improve the statistical and biological soundness of this analysis: switching from Brownian motion to an Ornstein-Uhlenbeck model for the phylogenetic transformation, and replacing Mahalanobis distance with PCA-based Euclidean distance.

---

## Improvement 1: Ornstein-Uhlenbeck phylogenetic correction

### Motivation

The current phylogenetic transformation assumes Brownian motion (BM), where trait variance increases linearly with divergence time. This is a neutral drift model. Physicochemical properties of functional proteins are not neutrally drifting — they are under stabilizing selection. A protein's molecular weight or isoelectric point is constrained by functional requirements, so distantly related orthologs don't diverge without bound.

Under BM, the phylogenetic VCV element between tips i and j is simply the shared evolutionary time:

```
C_BM[i,j] = t_MRCA(i,j)
```

Under the Ornstein-Uhlenbeck (OU) model, correlation decays exponentially with time since divergence:

```
C_OU[i,j] = (sigma^2 / 2*alpha) * exp(-alpha * (t_i + t_j - 2*t_MRCA)) * (1 - exp(-2*alpha * t_MRCA))
```

where alpha controls the strength of stabilizing selection. When alpha = 0, OU collapses to BM. When alpha is large, only very recently diverged proteins remain correlated.

**Practical consequence:** Under BM, the phylogenetic correction assumes that proteins separated by 500 million years should have vastly different physicochemical properties. When they don't (because of functional constraint), BM *over-corrects* — it attributes too much of the observed similarity to phylogeny, deflating the residual distances. Under OU, the model correctly expects that distant proteins may still be similar due to shared functional constraint, producing more accurate residual distances.

### Implementation

**Step 1: Estimate alpha for each gene family.**

Use `geiger::fitContinuous()` on one or more traits to estimate a single shared alpha per gene family. This adds roughly one line of R per gene family:

```r
# Fit OU model to estimate alpha (using first trait as representative)
ou_fit <- tryCatch(
  geiger::fitContinuous(gf_tree, gf_stats[, 1], model = "OU"),
  error = function(e) NULL
)
alpha <- if (!is.null(ou_fit)) ou_fit$opt$alpha else 0  # fallback to BM
```

Starting with a single shared alpha across all traits is the simplest approach and avoids a separate VCV per trait. Per-trait alpha or full multivariate OU (via `mvMORPH::mvOU()`) could be explored later but adds significant complexity.

**Step 2: Compute the OU VCV matrix.**

```r
ou_vcv <- function(tree, alpha) {
  C_bm <- ape::vcv.phylo(tree)
  t_root <- diag(C_bm)
  outer_times <- outer(t_root, t_root, "+")
  ou_C <- exp(-alpha * (outer_times - 2 * C_bm)) *
          (1 - exp(-2 * alpha * C_bm)) / (2 * alpha)
  dimnames(ou_C) <- dimnames(C_bm)
  return(ou_C)
}
```

**Step 3: Substitute into existing pipeline.**

In `phylo_gls_transform()` (`phylo_multivariate_distance_functions.R`), replace:

```r
phylo_vcv <- phytools::vcvPhylo(tree, anc.nodes = include_internal)
```

with:

```r
phylo_vcv <- ou_vcv(tree, alpha)
```

Everything downstream — the Cholesky decomposition in `phylo_correction.cpp`, the Mahalanobis computation — remains unchanged. The C++ code is agnostic to the source of the VCV matrix.

**Step 4: Handle edge cases.**

- If `fitContinuous` fails (common for very small gene families or trees with zero-length branches), fall back to BM (alpha = 0).
- If alpha is extremely large (> 100), cap it — this indicates the trait has essentially no phylogenetic signal, and the correction becomes numerically unstable.

### Files affected

- `bin/phylo_dist/phylo_multivariate_distance_functions.R` — add `alpha` parameter to `phylo_gls_transform()`, compute OU VCV
- `bin/phylo_dist/protein_distance_calculation_functions.R` — fit alpha in `calc_prot_spp_dists()`, pass to transform function
- No changes to C++ code (`phylo_correction.cpp`)

---

## Improvement 2: PCA-based Euclidean distance

### Motivation

After the phylogenetic transformation, the current pipeline computes pairwise Mahalanobis distances. This requires inverting the sample covariance matrix of the transformed data (a k x k matrix, where k = 10 physicochemical properties). Two problems arise:

**Problem 1: Small gene families.** Inverting a 10x10 covariance matrix requires at least 11 observations for the matrix to be full-rank. With `min_num_seq_per_og = 4`, many gene families have far fewer proteins than properties. The current code adds a tiny ridge regularization (`1e-6 * I` in `compute_inverse_covariance.cpp`), which prevents literal singularity but does not meaningfully stabilize the estimate. The resulting distances for small families are dominated by estimation noise.

**Problem 2: Correlated properties.** The 10 physicochemical properties are derived from the same amino acid sequence and are inherently correlated (e.g., molecular weight correlates with charge, GRAVY correlates with helix/sheet fractions). Mahalanobis distance accounts for this by standardizing through the covariance matrix, but with noisy covariance estimates, the standardization can amplify errors.

**The PCA solution:** Applying PCA to the transformed data and computing Euclidean distance on the resulting PC scores addresses both problems:

- PCA extracts at most `min(n-1, k)` meaningful components from n observations, naturally adapting to sample size. A gene family with 4 proteins yields 3 PCs, avoiding rank-deficiency entirely.
- By retaining only PCs explaining a threshold of total variance (e.g., 95%), noisy low-variance dimensions are discarded. This acts as principled regularization.
- The retained PCs are orthogonal (uncorrelated), so Euclidean distance is well-behaved — no covariance inversion needed.

**Relationship to Mahalanobis:** Euclidean distance on all PC scores, each scaled by `1/sqrt(eigenvalue)`, is mathematically identical to Mahalanobis distance. Unscaled Euclidean distance on PCs gives more weight to high-variance directions. After phylogenetic correction, the high-variance PCs capture the axes along which proteins diverge most beyond what phylogeny explains — arguably the most biologically relevant signal. Low-variance PCs are more likely estimation noise. Unscaled Euclidean on top PCs may therefore be both more robust and more biologically informative.

### Implementation

**Step 1: Replace `pairwise_mahalanobis()` with PCA + Euclidean distance.**

In `protein_distance_calculation_functions.R`, replace:

```r
dist_mat <- pairwise_mahalanobis(transf_data)
```

with:

```r
pca_result <- prcomp(transf_data, center = TRUE, scale. = FALSE)
cumvar <- cumsum(pca_result$sdev^2) / sum(pca_result$sdev^2)
n_pcs <- max(1, which(cumvar >= 0.95)[1])
scores <- pca_result$x[, 1:n_pcs, drop = FALSE]
dist_mat <- as.matrix(dist(scores))
```

**Step 2: Optionally output PC loadings for interpretability.**

PCA provides loadings that tell you which physicochemical properties drive the distance. This is an interpretability bonus over Mahalanobis, where the contribution of individual properties is opaque. Saving the loadings per gene family could help users understand *why* proteins are distant:

```r
loadings <- pca_result$rotation[, 1:n_pcs, drop = FALSE]
```

**Step 3: Remove or simplify C++ dependencies.**

The chunked Mahalanobis computation (`compute_inverse_covariance.cpp`, `compute_chunk_distances.cpp`, and the `process_all_chunks`/`get_dist_pair_chunks` helper functions in R) was necessary because Mahalanobis is computationally expensive for large matrices. Euclidean distance on a reduced-dimension matrix is trivial and can be computed with base R's `dist()`. These C++ files and the associated R helper functions can be removed or retained as dead code.

### Files affected

- `bin/phylo_dist/protein_distance_calculation_functions.R` — replace `pairwise_mahalanobis()` call with PCA + `dist()`
- `bin/phylo_dist/phylo_multivariate_distance_functions.R` — `pairwise_mahalanobis()`, `get_dist_pair_chunks()`, `process_single_chunk()`, `process_all_chunks()` can be removed
- `bin/phylo_dist/compute_inverse_covariance.cpp` — can be removed
- `bin/phylo_dist/compute_chunk_distances.cpp` — can be removed
- `modules/local/phylo_dist.nf` — no changes needed (output structure unchanged)

---

## Implementation order

The two improvements are independent and can be implemented in either order.

**Recommended order:**

1. **PCA + Euclidean first.** This is a simpler change (fewer lines of code, removes complexity rather than adding it) and addresses the most concrete statistical problem (covariance estimation failure in small gene families). It also simplifies the codebase by removing C++ dependencies.

2. **OU model second.** This is a more principled but less urgent change. BM over-correction is a systematic bias rather than a failure mode, and the PCA change partially mitigates its effects (by focusing on high-variance directions that are less affected by over-correction).

## Future considerations

- **Variance threshold as a parameter.** The 95% variance threshold for PC retention could be exposed as a pipeline parameter, allowing users to tune the regularization strength.
- **Per-trait alpha.** If certain physicochemical properties are under much stronger selection than others, per-trait alpha estimation would improve the OU correction. This requires a separate VCV per trait and a small refactor of the Cholesky transformation loop.
- **Multivariate OU.** For correlated traits under correlated selection, `mvMORPH::mvOU()` fits a full alpha matrix. This is the most principled model but computationally expensive per gene family. Worth investigating if results appear sensitive to the single-alpha assumption.
- **PC loadings for biological interpretation.** Outputting loadings per gene family enables downstream analyses asking *which* physicochemical properties drive protein divergence in each family.
