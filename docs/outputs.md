## Outputs

Note that a detailed walkthrough of how the results of NovelTree may be summarized and visualized, as applied to a dataset of 36 species of Telenemids, Stramenopiles, Alveolates, and Rhizarians, [may be found here](https://github.com/Arcadia-Science/2023-tsar-noveltree/tree/main/scripts/noveltree-summarization).

> **Mode-specific outputs**: Some outputs are only produced in certain workflow modes. These are marked with *(full mode only)* or *(zoogle mode only)*. Outputs without annotations are produced in all modes.

**1.** `preprocessing/` *(only with `--preprocess`)*: Preprocessed input proteomes, produced only when the pipeline is run with `--preprocess` (default off).

- `preprocessed_proteomes/`: One FASTA per species, named `{species}_preprocessed.fasta`. These are the full-length, unmasked protein sequences after quality cleanup (trailing stop codons removed; internal stops and rare residues recoded), minimum-length filtering (`--min_protein_length`), and — for non-reference proteomes — CD-HIT deduplication (identity 1.00, or 0.97 for TransDecoder-derived proteomes). Reference proteomes skip CD-HIT and retain all sequences. Written via Nextflow `storeDir`, so this directory persists and is reused as a cache across runs.

**2.** `busco/` *(full mode only)*: Contains output of all BUSCO analyses. Primary directory includes:  

- Short summary results for each species, to each lineage dataset (shallow or broad) in text and json formats.  
- Batch summary results for each species, to each lineage dataset.  
- Directory for each species' busco analysis with more detailed output.  

**3.** `protein_annotations/`: InterPro domain annotations obtained per each species for which UniProt accessions are available. Used by COGEQC to assess orthogroup quality during MCL inflation parameter selection.

**4.** `diamond/`: contains the results of all pairwise comparisons of sequence similarity, between and within species using diamond BLASTP.  

**5.** `orthofinder/`: contains all results from OrthoFinder runs, including the MCL testing stage (`mcl_test_dataset`) and the full analysis using the best-performing inflation parameter (`complete_dataset`).

- `mcl_test_dataset/`: Contains one directory per tested inflation parameter. Directory structure follows OrthoFinder convention - gene family sequences are not retained to save space.
- `complete_dataset/`: OrthoFinder output for the best-performing inflation parameter following standard convention (using OrthoFinder's -os flag); holds the `Results_Inflation_<value>/` directory.
- Published directly under `orthofinder/` (alongside `complete_dataset/`):
  - `all_ogs_counts.csv`: comma-separated CSV listing, for all orthogroups, the number of included species, total copy number, mean copy number, and number of higher-level taxonomic groups included.
  - `spptree_core_ogs_counts.csv` / `genetree_core_ogs_counts.csv`: the same as above, but for only the two respective subsets of gene families (species-tree set: high species coverage, low copy number; gene-tree set: broader, ≥4 species). These CSVs define which gene families enter each downstream analysis.
  - `chimera_report.tsv`: report of cross-orthogroup chimeric proteins flagged and removed.

  The per-orthogroup *unaligned* FASTAs are published flat under `alignments/unaligned/` (see section 7), not here.

- `mcl/`: Nextflow `storeDir` resume cache for the MCL step. It mirrors the same `complete_dataset/` and `mcl_test_dataset/` directories (and the CSVs above) that are published directly under `orthofinder/`; these copies exist to enable `-resume` and can be safely deleted to reclaim space.

**6.** `orthogroup_summaries/`: Results from COGEQC.

- A table reporting all calculated summary statistics for the orthogroups inferred for each inflation parameter.
- A text file that lists the best-performing inflation parameter.
- A PDF plotting the summary statistics as a function of inflation parameter value.

**7.** `alignments/`: Per-gene-family sequences organized by processing stage (unaligned → aligned → masked → trimmed).

- `unaligned/`: Flat directory of the unaligned OrthoFinder per-orthogroup FASTAs (one per gene family) that enter the alignment step. Limited to the families used in phylogenetic analysis (the union of the species-tree and gene-tree core sets); small/taxon-specific families excluded from analysis are not included.
- `original/`: Aligned-but-unmasked output from the selected aligner, named by aligner: MAFFT tier-1 writes `{OG}_einsi.fa` (or `{OG}_linsi.fa` when `--mafft_mode linsi`; the same mode is applied to all tier-1 families), FAMSA writes `{OG}_famsa.fa`, and the legacy single-aligner MAFFT path writes `{OG}_mafft.fa`. For WITCH this is the pre-masking alignment, saved as `{OG}_witch_unmasked.fa`.
- `masked/` *(WITCH only)*: WITCH's confidence-masked and length/gap-filtered alignments (`{OG}_witch.fa`). This is the alignment that feeds downstream gene-family tree inference.
- `trimmed/`: Externally trimmed alignments from ClipKIT (`{OG}_*_clipkit.fa`) or CIAlign (`{OG}_*_cialign.fa`), produced only when an external trimmer is enabled.
- `species_protein_maps/` (within each stage subdirectory, e.g. `original/`, `masked/`, `trimmed/`): Protein-to-species mapping files (`*_map.link`) used by downstream reconciliation tools.

**8.** `gene_family_trees/`: Gene family trees organized by processing stage.

- `original/`: Gene family trees as inferred by FastTree2 or IQ-TREE. Files are named `{alignment_prefix}_{method}.newick` where method is `ft` (FastTree) or `iqt` (IQ-TREE), and the alignment prefix encodes the aligner and trimmer used (e.g., `OG0000001_witch_clipkit_iqt.newick`).
- `reconciled/`: Gene family trees reconciled with the species tree by GeneRax. Contains `generax_per_species/` (per-species rate model) and `generax_per_family/` *(full mode only)* (per-family rate model). Files are named `{OG}_reconciled_gft.newick`.
- `time_calibrated/` *(zoogle mode only)*: Time-calibrated gene family trees dated using speciation node ages from the species tree.
  - `{OG}_dated.newick`: Time-calibrated gene family tree.
  - `{OG}_calibrations.csv`: Speciation-node calibrations used for dating, with columns: `gf_mrca` (gene tree node), `age_mya` (species tree age), `tipA`/`tipB` (MRCA-defining tips), `n_desc_tips`, `events_node`, `min_mya`, `max_mya` (±`age_bracket`, default ±20%).

**9.** `species_trees/`: Species trees organized by inference method.

- `asteroid/`: Unrooted (or outgroup-rooted) species tree from [Asteroid](https://github.com/BenoitMorel/Asteroid). Includes `asteroid.bestTree.newick` (best-scoring tree), `asteroid.allTrees.newick`, `asteroid.bsTrees.newick` (bootstrap trees), `asteroid.scores.txt`, and `disco_decomposed_rooted_gfts.newick` (single-copy trees from [DISCO](https://github.com/JSdoubleL/DISCO) decomposition).
- `speciesrax/`: Rooted species tree from [SpeciesRax](https://github.com/BenoitMorel/GeneRax/wiki/SpeciesRax). Contains inferred species trees with support values (`inferred_species_tree.newick`, `species_tree_*.newick`, `starting_species_tree.newick`), per-species event and coverage summaries, and run statistics.
- `time_calibrated/` *(zoogle mode only)*: Time-calibrated species tree. Contains `time_calibrated_species_tree.newick` and `calibration_log.txt`.

**10.** `generax/`: Results from GeneRax gene-family tree/species tree reconciliation with joint optimization of gene tree topology and duplication/loss rates. Each gene family directory contains the key output files in a flat structure, plus a compressed archive of the full GeneRax output.

- `per_family_rates/` *(full mode only)*: Results using the per-family rate model. One directory per gene family containing:
  - `{OG}_reconciled_gft.newick`: Reconciled gene family tree.
  - `{OG}_eventCounts.txt`: Counts of duplication and loss events.
  - `{OG}_speciesEventCounts.txt`: Event counts per species-tree node.
  - `{OG}_full_output.tar.gz`: Compressed archive of the complete GeneRax output for this gene family.
- `per_species_rates/`: Results using the per-species rate model. One directory per gene family containing the same files as above, plus:
  - `{OG}_perSpeciesCoverage.txt`: Per-species coverage statistics.

**11.** `gene_family_evolution/`: Phylogenetic profiles and event summaries from PHYLO_PROFILES module.

- `duplication_count_per_species_per_gene_family.tsv`: Matrix of duplication events with species-tree nodes as rows and gene families as columns.
- `loss_count_per_per_species_gene_family.tsv`: Matrix of gene loss events per species-tree node per gene family.
- `speciation_count_per_species_per_gene_family.tsv`: Matrix of speciation events per species-tree node per gene family.

**12.** `protein_properties/` *(zoogle mode only)*: Amino acid composition and physicochemical property summaries from PROTEIN_PROPERTIES module.

- `aa-summary-stats/across-family-summaries/`: Aggregated statistics across all gene families.
- `aa-summary-stats/per-family-summaries/`: Per-family directories containing property summaries for each protein.
  - Contains CSV files with 20 amino acid frequencies and physicochemical properties (molecular weight, aromaticity, instability, flexibility, GRAVY, isoelectric point, charge at pH 3/5/7/9, helix/sheet fractions, extinction coefficients).

**13.** `zoogle/` *(zoogle mode only)*: Phylogenetically-corrected protein distance analyses from ZOOGLE module. Produced for all gene families with at least 4 proteins from at least 2 species (each with ≥2 proteins).

*Universal outputs (always produced):*

- `phylo-corrected-data/`: Physicochemical property data after phylogenetic correction.
- `protein-dist-mats/`: Raw protein-protein distance matrices based on physicochemical properties.
- `protein-phylo-dist-mats/`: Phylogenetically-corrected Mahalanobis distance matrices.
- `centroid-dists/`: Mahalanobis distance from each protein to the family centroid (column means of GLS-transformed data). Includes rank and empirical protein-level p-values.
- `centroid-summary-tables/`: Comprehensive centroid summary tables with gene family, protein, species, centroid distance, rank, protein-level p-value, and species-level permutation p-value.

*Reference-specific outputs (only when `ref_species` is set and present in the gene family):*

- `protein-dists-to-reference/`: Distances from each protein to reference species proteins.
- `species-dists-to-reference/`: Species-level average distances to reference.
- `protein-pvals/`: Statistical significance (p-values) of protein distances to reference species.
- `species-pvals/`: Statistical significance (p-values) of species-level distances to reference.
- `pairwise-protein-dist-perm-test/`: Permutation test results for pairwise protein comparisons.
- `final_protein_pair_summary_tables/`: Comprehensive summary tables combining distance metrics, p-values, and protein metadata for each gene family.

**14.** `pipeline_info/`: Nextflow execution reports and pipeline metadata.

- `execution_timeline_*.html`: Interactive timeline showing when each process started and completed.
- `execution_report_*.html`: Detailed execution report with resource usage statistics.
- `execution_trace_*.txt`: Tab-separated file with detailed metrics for each task.
- `pipeline_dag_*.html`: Directed acyclic graph visualization of the workflow.  
