# Arcadia-Science/noveltree: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0-alpha - 2026-03-30

### Added
- Zoogle analysis mode (`-profile zoogle`): end-to-end physicochemical protein distance analysis with time-calibrated gene family trees, phylogenetic correction, Mahalanobis distances, and permutation tests
  - PROTEIN_PROPERTIES, DATE_GENE_FAMILY_TREES, TIME_CALIBRATE_SPECIES_TREE, BUILD_REFERENCE_CHRONOGRAM, ZOOGLE_ANALYSIS modules
  - Centroid-based analysis mode (`--ref_species none`) for analysis without a reference species
- Simplified execution mode (`-profile simplified`): streamlined pipeline for large datasets (adaptive alignment default, no BUSCO, per-species GeneRax EVAL strategy only)
- Arcadia production profile (`-profile arcadia`): combines zoogle + AWS Batch
- Adaptive three-tier alignment routing (`--aligner adaptive`): MAFFT_TIER1 (≤200 seqs) → WITCH_TIER2 (≤3000) → FAMSA_TIER3 (larger), with automatic fallback between tiers
- FAMSA aligner option (`--aligner famsa`)
- IQ-TREE with automatic FastTree fallback (`--iqtree_fasttree_fallback`)
- PARSE_PHYLOHOGS module: per-OG ortholog/paralog inference and HOG membership directly from GeneRax NHX reconciliation (`extract_relationships_from_nhx.py`)
- PHYLO_PROFILES module for summarizing gene duplication, loss, and speciation events per species per gene family
- Optional proteome preprocessing (`--preprocess`): TransDecoder for transcriptomes, isoform filtering, minimum protein length filtering
- Samplesheet support for UniProt proteome IDs (`UP*`) and NCBI genome accessions (`GCA_*/GCF_*`) as data sources, with automatic download
- `--test_run_mcl` flag for optional MCL inflation parameter testing on species with UniProt accessions
- Dynamic resource allocation for GeneRax, WITCH, IQ-TREE, and MAFFT based on gene family size (sequence count and max length)
- RENAME_FASTAS process to standardize FASTA filenames before all downstream processes
- AWS Batch support (`-profile awsbatch`) with `--awsqueue` and `--awsregion` parameters
- Singularity support (`-profile singularity`) with automatic Docker-to-Singularity conversion and image caching
- `--test_run` flag: restrict analysis to gene families containing all species (fast smoke test)
- Eukaryote test dataset (6 Opisthokont species)
- `Makefile` for Docker image builds, `CITATIONS.md`, `docs/singularity.md`
- New dependencies: treePL, PATHd8, phangorn

### Changed
- Main workflow refactored from monolithic `main.nf` (~500 lines) to 7 focused subworkflows (INPUT_CHECK, PREPARE_INPUTS, INFER_ORTHOGROUPS, INFER_GENE_TREES, RECONCILE_TREES, RECONCILIATION_SUMMARIES, ZOOGLE)
- Default aligner: `witch` → `adaptive` (three-tier routing enabled by default)
- Default tree method: `fasttree` → `iqtree`
- Default MSA trimmer: `none` → `clipkit`
- `min_ungapped_length` default: 20 → 50
- `min_num_spp_per_og` default: 4 → 2
- `min_prop_spp_for_spptree` default: 0.25 → 0.75
- `max_copy_num_spp_tree` default: 5 → 10
- `mcl_inflation` default: `'1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0'` → `'1.5'`
- FastTree arguments updated to evidence-based maximum accuracy settings (`-lg -gamma -bionj -pseudo -spr 4 -mlacc 2 -slownni`; Zhou et al. 2018)
- GeneRax container updated to `generax_56f3ed0:1.1.3`
- GeneRax reconciliation model: `UndatedDTL` → `UndatedDL` (gene transfers removed)
- Rare amino acid handling: selenocysteine (U) and pyrrolysine (O) recoded as X before alignment and reconciliation
- Species name handling standardized: hyphens used throughout (underscores/spaces auto-converted)
- UniProt annotation retrieval rewritten to use ID Mapping API (replaces bioservices)
- MCL inflation selection uses only InterPro scoring (OMA removed — Sorensen-Dice incompatible with 1:1 OMA group IDs)
- Gene family tree dating uses speciation-only calibrations from GeneRax reconciliation with treePL fixed smoothing (smooth=10) and `--age_bracket` (default ±20%) calibration brackets (replaces congruification)
- FILTER_ORTHOGROUPS functionality absorbed into ORTHOFINDER_MCL
- `MAFFT_ADAPTIVE` renamed to `MAFFT_TIER1` for consistency with WITCH_TIER2 and FAMSA_TIER3
- Samplesheet redesigned: required columns reduced from 7 to 3 (`species`, `input_data`, `input_type`); all other columns are optional in any order. Column renames for clarity: `file` → `input_data`, `mode` → `input_type`, `uniprot` → `has_uniprot_ids`, `mcl_test` → `include_in_mcl_test`, `isoform` → `filter_isoforms`, `reference` → `reference_proteome`, `shallow_db` → `busco_shallow`, `broad_db` → `busco_broad`
- All samplesheet boolean columns standardized to `yes`/`no` (previously `uniprot` used `true`/`false`)

### Removed
- ORTHOFINDER_PHYLOHOGS module (replaced by PARSE_PHYLOHOGS using GeneRax NHX directly)
- FILTER_ORTHOGROUPS module (absorbed into ORTHOFINDER_MCL)
- `bin/translate_gene_trees.py` (was used by ORTHOFINDER_PHYLOHOGS)
- PMSF two-pass tree inference (`iqtree_pmsf.nf`)
- `min_num_grp_per_og`, `max_copy_num_gene_trees`, `tree_model_pmsf` parameters
- `species_tree_prep` module
- OMA annotation collection and scoring from cogeqc analysis
- `bin/protein_annotation.R` (dead code; only Python version was used)
- `taxonomy` samplesheet column and taxonomy group counting from orthogroup summaries
- `bioservices` dependency (replaced by `requests` for UniProt ID Mapping API)
- Congruification approach for gene family tree dating (replaced by reconciliation-filtered speciation-only calibrations)
- Gene transfer (HGT) inference from GeneRax reconciliation (UndatedDTL → UndatedDL)

## v1.0.1-alpha - 09/28/2023

Release of NovelTree that is associated with the pub ["NovelTree: Highly parallelized phylogenomic inference"](https://doi.org/10.57844/arcadia-z08x-v798). Includes small bug fix caused by including the `xref_tigrfam` return field when querying UniProt. 

## v1.0.0-alpha - 09/27/2023

Initial release of NovelTree. Do not use - this version is deprecated in favor of v1.0.1-alpha.
