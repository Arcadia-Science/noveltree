# Arcadia-Science/noveltree: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0-alpha - 2026-03-06

### Added
- Zoogle analysis mode (`-profile zoogle`): end-to-end physicochemical protein distance analysis with time-calibrated gene family trees, phylogenetic correction, Mahalanobis distances, and permutation tests
  - ZOOGLE, PHYSICOCHEMICAL_PROPS, DATE_GENE_FAMILY_TREES, TIME_CALIBRATE_SPECIES_TREE, BUILD_REFERENCE_CHRONOGRAM modules
- Simplified execution mode (`-profile simplified`): streamlined pipeline for large datasets (FAMSA, no BUSCO, per-species GeneRax only)
- FAMSA aligner option (`--aligner famsa`)
- IQ-TREE with FastTree fallback (`--iqtree_fasttree_fallback`)
- PHYLO_PROFILES and MERGE_PHYLO_PROFILES modules for summarizing gene duplication, transfer, loss events
- ORTHOFINDER_PHYLOHOGS module for hierarchical ortholog inference
- RENAME_FASTAS process to standardize FASTA filenames before all downstream processes
- AWS Batch support (`-profile awsbatch`) with `--awsqueue` and `--awsregion` parameters
- Singularity support (`-profile singularity`) with automatic Docker-to-Singularity conversion
- Eukaryote test dataset (6 Opisthokont species, 50 orthogroups)
- `Makefile` for Docker image builds, `CITATIONS.md`, `docs/singularity.md`

### Changed
- Default MSA trimmer: `none` → `clipkit`
- `max_copy_num_spp_tree` default: 5 → 10
- Aligner and tree inference refactored into subworkflows
- GeneRax container updated to `generax_56f3ed0:1.1.3`
- Improved rare amino acid handling (U→X, O→X before alignment)
- Standardized species name handling: hyphens used throughout (underscores/spaces auto-converted)
- UniProt annotation retrieval rewritten to use ID Mapping API (replaces bioservices)
- MCL inflation selection uses only InterPro scoring (OMA removed — Sorensen-Dice incompatible with 1:1 OMA group IDs)
- Gene family tree dating uses speciation-only calibrations from GeneRax reconciliation (replaces congruification)

### Removed
- PMSF two-pass tree inference
- `min_num_grp_per_og`, `max_copy_num_gene_trees`, `tree_model_pmsf` parameters
- `species_tree_prep` module
- OMA annotation collection and scoring from cogeqc analysis
- `bin/protein_annotation.R` (dead code; only Python version was used)
- `bioservices` dependency (replaced by `requests` for UniProt ID Mapping API)
- Congruification approach for gene family tree dating

## v1.0.1-alpha - 09/28/2023

Release of NovelTree that is associated with the pub ["NovelTree: Highly parallelized phylogenomic inference"](https://doi.org/10.57844/arcadia-z08x-v798). Includes small bug fix caused by including the `xref_tigrfam` return field when querying UniProt. 

## v1.0.0-alpha - 09/27/2023

Initial release of NovelTree. Do not use - this version is deprecated in favor of v1.0.1-alpha.
