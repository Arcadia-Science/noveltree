## Outputs
Note that a detailed walkthrough of how the results of NovelTree may be summarized and visualized, as applied to a dataset of 36 species of Telenemids, Stramenopiles, Alveolates, and Rhizarians, [may be found here](https://github.com/Arcadia-Science/2023-tsar-noveltree/tree/main/scripts/noveltree-summarization).  

**1.** `busco/`: Contains output of all BUSCO analyses. Primary directory includes:  

- Short summary results for each species, to each lineage dataset (shallow or broad) in text and json formats.  
- Batch summary results for each species, to each lineage dataset.  
- Directory for each species' busco analysis with more detailed output.  

**2.** `protein_annotations/`: protein annotations obtained per each species for which uniprot accessions (i.e. those corresponding to RefSeq protein accessions) are available.  

- Which annotations are included here depends on how the "download_annots" parameter was specified.  

**3.** `diamond/`: contains the results of all pairwise comparisons of sequence similarity, between and within species using diamond BLASTP.  

**4.** `orthofinder/`: contains all results from orthofinder runs, including the MCL testing stage (`mcl_test_dataset`) and the full analysis using the best-performing inflation parameter (complete_dataset).  

- `mcl_test_dataset/`: Contains one directory per tested inflation parameter. Directory structure follows OrthoFinder convention - gene family sequences are not retained to save space.  
- `complete_dataset/`: Contains OrthoFinder output for the best-performing inflation parameter following standard convention (using OrthoFinder's -os flag).  

**5.** `orthogroup_summaries/`: Results from COGEQC.  

- A table reporting all calculated summary statistics for the orthogroups inferred for each inflation parameter.  
- A text file that lists the best-performing inflation parameter.  
- A PDF plotting the summary statistics as a function of inflation parameter value.  

**6.** `filtered_orthogroups/`: final orthogroups used for gene-family tree and species tree inference that passed user-specifed filters.  

- `species_tree_og_msas/`: trimmed multiple sequence alignments for gene families used in species tree inference.  
- `gene_tree_og_msas/`: trimmed multiple sequence alignments for gene families used only in gene-family tree / species tree reconciliation.  
- `all_ogs_counts.csv`: comma-separated csv listing, for all orthogroups (including those fro which msa/gene family trees are not inferred), the number of included species, total copy number, mean copy number, and number of higher-level taxonomic groups included.  
- `(gene)speciestree_core_ogs_counts.csv`: the same as above, but for only the two respects subsets of gene families.

**7.** either `mafft_alignments/` or `witch_alignments`: multiple sequence alignments for orthogroups passing filtering thresholds (e.g. minimum number of species, maximum copy number).   
- `witch_alignments` contains both the "raw" alignments inferred from witch (`original_alignments`), as well as the cleaned alignments produced by the software, with gappy or otherwise poorly alignmed columns removed.  
- If no alignment cleaning method is used, these directories will contain another subdirectory, `species_protein_maps`, which contains the map files linking each protein ID to the parent species.  

**8.** `trimmed_msas/`: Only produced if using CIAlign or ClipKIT to clean multiple sequence alignments.  
- If this directory is produced by the workflow, the `species_protein_maps` subdirectory can be found here.   

**9.** `fasttree_gene_trees` or `iqtree_gene_trees/`: Gene family trees using either FastTree2 or IQ-TREE respectively.  

**10.** `asteroid/`: Asteroid species tree, either rooted or unrooted. Includes:  

- asteroid.bestTree.newick: Single species tree with the greatest likelihood among set of inferred trees (for instance if multiple random starting trees are used).  
- asteroid.allTrees.newick: All inferred species trees (=1 if default of 1 random starting tree is used).  
- asteroid.bsTrees.newick: Bootstrapped species trees. Number of trees = number of bootstrap replicates.   
- asteroid.scores.txt: likelihood scores for all inferred trees.  
- disco_decomposed_rooted_gfts.newick: A newick tree file containing all single-copy gene family trees inferred using [DISCO](https://github.com/JSdoubleL/DISCO) by decomposing each mutiple copy gene family tree into their respective single-copy counterparts.  

**12.** `speciesrax/`: all results/outputs from SpeciesRax inferred using the subset of gene families that passed filters for involvement in rooted species tree inference. Full description of these outputs (including gene-family tree/species tree reconciliations) are described on the GeneRax github. Key output directories includes:  

- `species_trees/`: contains inferred rooted species trees, species tree likihoods, and species trees with internal nodels labeled according to their support values.  
- `reconciliations/`: gene-family tree - species tree reconciliations (species trees with gene family duplications, transfers, and losses mapped on) in several formats. May be plotted using reconciliation software like thirdkind. Additionally contains files for each gene family that enumerate duplication-transfer-loss event counts per species. Reconciliations here are obtained without optimization of the gene family tree topology under a model of gene duplication, transfer, and loss.  

**13.** `generax/`: all results/outputs from GeneRax inferred using the subset of gene families that passed filters for involvement in rooted species tree inference. Full description of these outputs (including gene-family tree/species tree reconciliations) are described on the GeneRax github. Reconciliations are the result of joint optimization of the gene family tree topology and reconciliation, and duplication/transfer/loss rates.  

- `per_family_rates`: Results of the per-family model. One directory for each gene family, named by the corresponding family ID.  
- `per_species_rates`: Results of the per-species model. One directory for each gene family, named by the corresponding family ID.  
- Each has a similar directory structure to SpeciesRax, but lacking the species trees directory, and lacking another directory containing the results of gene family tree reconciation and inferred rates of gene duplication, transfer, and loss.  
- `results/`: one directory per-gene-family containing reconciled gene trees and inferred rates of gene family duplication, transfer and loss.  
