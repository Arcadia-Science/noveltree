# Run the workflow on your own custom dataset
As a reminder - the pipeline assumes that prior to analysis, input proteomes for each species have been sufficiently filtered such that no additional filtering of species or sequences is required. For a description of such a filtering procedure, see the following [GitHub repository](https://github.com/Arcadia-Science/2023-tsar-noveltree/tree/main/scripts/data-preprocessing).

## Preparation
**1.** Prepare a samplesheet following the required format:  

```
species,file,taxonomy,shallow_db,broad_db,mode,uniprot,mcl_test
Entamoeba_histolytica,Entamoeba_histolytica-test-proteome.fasta,Amoebozoa,NA,eukaryota_odb10,proteins,true,true
```

> #### Description of Columns:
>
> `species`: species name to use.  
> `file`: complete path to fasta file, whether local or remote (e.g. provide complete local file path, or S3 URI/hyperlink to other cloud storage).  
> `taxonomy`: higher-rank taxonomy for each species (e.g. supergroup, class, family, genus). Utility of this parameter depends on the taxonomic scope of each dataset. Used in filtering orthogroups for phylogenetic inference.  
> `shallow_db`: busco lineage dataset for shallow taxonomic scale analysis (e.g. below eukaryota). If NA, will not run.  
> `broad_db`: busco lineage dataset for broad taxonomic scale analysis (e.g. eukaryota). If NA, will not run.  
> `mode`: specification of busco analysis mode.  
> `uniprot`: true/false specification indicating whether the proteome comes from UniProt (i.e. has UniProt protein accessions that `NovelTree` can use to annotate).  
> `mcl_test`: true/false specification of whether this species is to be included in the MCL inflation parameter test-set. These species must have UniProt protein accessions (for COGEQC protein domain score).  

**2.** Create a parameter file that includes all necessary input, output, and parameter specifications: example below.  

```
{
  "input": "/full/path/to/samplesheet.csv",
  "outdir": "results",
  "mcl_inflation": "1.0,2.0,3.0",
  "min_ungapped_length": 30,
  "min_num_spp_per_og": 4,
  "min_num_grp_per_og": 1,  
  "aligner": "witch",
  "max_copy_num_spp_tree": 5,
  "max_copy_num_gene_trees": 10,
  "download_annots": "none",
  "tree_model": "LG+F+G4",
  "outgroups": "none"
}
```

> #### Parameter descriptions:
>
> `input`: Complete filepath to input samplesheet. May be locally stored, or remotely stored (again - if remote, provide S3 URI, or hyperlink to other cloud storage).  
> `mcl_inflation`: DEFAULT "1.5,2.0,2.5,3.0". Quoted, comma-separated list of MCL inflation parameters to be tested when clustering proteins into orthogroups with OrthoFinder.  
> `min_ungapped_length`: DEFAULT: 20. The minimum ungapped length of cleaned/trimmed multiple sequence alignments.  
> `min_num_spp_per_og`: DEFAULT: 4. Minimum # of species a gene family must contain for phylogenetic inference.  
> `min_num_grp_per_og`: DEFAULT: 1. Minimum # of 'higher' order taxonomic groups an gene family must contain for phylogenetic inference.  
> `aligner`: DEFAULT: "witch". Method used to infer multiple sequence alignments. Either MAGUS ('magus') or MAFFT ('mafft').  
> `max_copy_num_spp_tree`: DEFAULT: 5. Maximum # of per-species gene copy number a gene family may contain for species-tree inference.  
> `max_copy_num_gene_trees`: DEFAULT: 10. Maximum # of per-species gene copy number a gene family may contain for gene tree - species tree reconciliation with GeneRax.  
> `min_prop_spp_for_spptree`: DEFAULT: 0.25. Minimum proportion of species a gene family must contain to be used in species tree inference.  
> `download_annots`: DEFAULT: "minimal". Set of annotations to be downloaded. "none" corresponds to a minimal set. See description of parameters for expanded description of options.  
> `tree_model`: DEFAULT: "LG+F+G4". Model of amino acid substition to be used for phylogenetic inference. If using a posterior mean site frequency model (see below), this model will be used to infer an initial guide-tree.  
> `tree_model_pmsf`: OPTIONAL: Posterior mean site frequency model to be used for phylogenetic inference (e.g. "LG+C40+F+G4"). If not specified (i.e. excluded from parameter file), only `tree_model` will be used.  
> `outgroups`: OPTIONAL: A comma separated string of species IDs to be used to manually root Asteroid species tree. If specified, this species tree will have branch lengths estimated with SpeciesRax, and will be used for all GeneRax analyses.  
>
> Alternatively, you can use the test dataset provided by Arcadia Science [here](https://github.com/Arcadia-Science/test-datasets/noveltree).  

**3.** Ensure that proteins are named following the following convention: `Species_genus:GeneID`  

```
# Example:
>Entamoeba_histolytica:C4M6M9 NGG1-interacting factor 3

# Everything prior to the colon (:) is a constant identifier unique to that species/lineage
# and can be whatever you would like, but must not include a colon.

# Everything that follows the colon is what must be a unique protein/gene identifier
# Additional sequence info may be included following a space.

# If you intend to download annotations for a proteome, the sequence identifier must
# be the UniProt protein accession. NovelTree uses the string that follows the colon
# to extract the uniprot accession and annotate proteins.

# Future versions will include a utility to automate sequence naming, and the ability
# to automatically correspond other standard sequence identifiers (e.g. NCBI RefSeq)
# with UniProt accessions to facilitate this annotation process.
```

**4.** Start running your own analysis! See [here](#parameter-specification) for in-depth parameter description.  

```bash
nextflow run . -profile docker -params-file <PARAMS.JSON>
```

---

### **NOTE: Regarding our analysis on Nextflow Tower**  

When NovelTree was applied to the dataset used within its associated pub, we did so on Nextflow Tower. Accordingly, we specified a handful of additional configurations.  

This included:  
  1) `max_cpus = 5000`: This set the maximum number of available cpus (as spot instances) to all concurrent processes. Effectively the number of CPUs available to our virtual "cloud" computer.  
  2) `max_memory = 30000.GB`: The same, but for memory alloted for all concurrent processes.  
  3) `max_time = 2400.h`: Again, the same, but the maximum time alloted for all concurrent processes.  
  4) Additionally, we allocated 32 CPUs to the head node to ensure efficient monitoring and submission of jobs.

![Workflow Figure](https://github.com/Arcadia-Science/noveltree/blob/ap/readmes/NovelTree.png)  

## The workflow proceeds to conduct the following steps:
1. `INPUT_CHECK`: Proteomes are staged locally (including downloaded from S3 or other cloud storage if necessary)  
2. `BUSCO`: Each proteome is summarized using [`BUSCO`](https://busco.ezlab.org/) completeness at both user-specified shallow (e.g. Eukaryota) and taxon-specific scales  
3. `PROTEIN_ANNOTATION`: Proteomes for which sequence names include [`UniProt`](https://www.uniprot.org/) protein accessions are annotated using [`UniProt.ws`](https://bioconductor.org/packages/release/bioc/html/UniProt.ws.html)  
4. `ORTHOFINDER_PREP`: Proteomes are staged/reformated for analysis with [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)  
5. `DIAMOND_BLASTP`: Determine all-v-all (within and among species) protein sequence similarity using [`Diamond`](https://github.com/bbuchfink/diamond) BlastP ultra-sensitive  
6. `ORTHOFINDER_MCL`: Cluster [`UniProt`](https://www.uniprot.org/) sequences into orthogroups/gene-families using [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)'s implementation of [`MCL`](http://micans.org/mcl/) clustering using a specified set of inflation scores  
7. `COGEQC`: Summarization and quantification of gene family inference performance using a set of summary statistics, including the functional annotation score using [`COGEQC`](https://almeidasilvaf.github.io/cogeqc/index.html) applied to both [`InterPro`](https://ebi.ac.uk/interpro/) domain annotations and [`OMA`](https://omabrowser.org/oma/home/) orthology IDs.  
8. `SELECT_INFLATION`: Based on the above summaries, select the (mean) inflation parameter that performs best (e.g. orthogroups are most homogenous in protein domain annotations, penalizing against dispersal of annotations across orthogroups), accounting for diminishing returns with increasing or decreasing parameter values.  
9. `ORTHOFINDER_MCL`: Repeat step six (6: MCL clustering into orthogroups) for all species under the optimal inflation parameter  
10. `FILTER_ORTHOGROUPS`: Summarize distribution of orthogroups across taxonomic groups and per-species copy number, filtering into a conservative subset for species tree inference, and one for gene-family tree inference.  
11. `ALIGN_SEQS`: Infer multiple sequence alignments for each focal gene family with either [`WITCH`](https://github.com/c5shen/WITCH) (default) or [`MAFFT`](https://mafft.cbrc.jp/alignment/software/)  
12. `TRIM_SEQS`: OPTIONAL: Trim uninformative/memory-consuming/gappy segments of alignments with either [`CIAlign`](https://github.com/KatyBrown/CIAlign) (defualt) or [`ClipKit`](https://jlsteenwyk.com/ClipKIT/)  
14. `INFER_TREES`: Infer gene family trees using either [`FastTree2`](http://www.iqtree.org/) (default) or [`IQ-TREE`](http://www.iqtree.org/)  
16. `ASTEROID`: Infer an unrooted species tree using [`Asteroid`](https://github.com/BenoitMorel/Asteroid). If outgroups are specified, this tree will be rooted using these species.  
17. `SPECIESRAX`: Infer a rooted species tree, estimating its topology under a model of gene duplication, transfer, and loss using [`SpeciesRax`](https://github.com/BenoitMorel/GeneRax/wiki/SpeciesRax). If outgroups are provided, [`SpeciesRax`] infers branch lengths for the `ASTEROID` tree.  
18. `GENERAX_PER_FAMILY`: Reconcile gene family trees with the species tree, inferring rates of gene duplication, transfer and loss using [`GeneRax`](https://github.com/BenoitMorel/GeneRax) under the per-family model (rates are constant across all species/branches)  
18. `GENERAX_PER_SPECIES`: Reconcile gene family trees with the species tree, inferring rates of gene duplication, transfer and loss using [`GeneRax`](https://github.com/BenoitMorel/GeneRax) under the per-species model (each species/branch has own rates)  
19. `ORTHOFINDER_PHYLOHOGS`: Infer phylogenetically hierarchical orthologs using [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)  
