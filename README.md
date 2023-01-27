## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**Arcadia-Science/phylorthology** is a phylogenomic pipeline designed to proteomes from diverse organisms and inferring orthology, gene-family trees, and a species tree.

The pipeline assumes that prior to analysis, input proteomes for each species have been sufficiently filtered such that no additional filtering of species or sequences is required. (NOTE: May still explore species-level filtering based on MCL clustering? e.g. for each species, count their # of single species orthogroups for each inflation parameter?)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.

## Pipeline summary
At its core, PhylOrthology used the OrthoFinder algorithm to normalize pairwise protein similarity scores to account for evolutionary divergence prior to clustering into orthogroups/gene families with MCL clustering. Because this clustering is contingent upon the MCL inflation parameter, PhylOrthology automates the identification of the inflation parameter that returns the most biologically sensible and tractable set of orthogroups.

Thus, two rounds of protein clustering takes place - an initial round for inflation parameter testing on a (reduced) set of proteomes for which UniProt protein accessions are available, and a second round on the complete dataset. Because Nextflow can automate the parallelization of tasks that are independent of each other, PhylOrthology will begin the all-v-all protein comparisons for both rounds simultaneously across available compute resources. Once the first round of MCL clustering has completed, we summarize orthogroups based on a number of metrics, choosing a best-performing inflation parameter for the analysis of the full dataset. This includes a functional protein domain annotation score calculated with COGEQC, which quantifies the ratio of InterPro domain "Homogeneity" of domains within orthogroups to "Dispersal" of domains among orthogroups.

With orthogroups/gene families inferred, we will summarize each orthogroup on the basis of their taxonomic and copy number distribution, quantifying the number of species/clades included in each, as well as the mean per-species copy number. These summaries facilitate 'filtering' for sufficiently conserved/computationally tractable gene families for downstream phylogenetic analysis. In other words, it's best to avoid excessively small (e.g. < 4 species) or large gene families (e.g. > 200 species and mean copy # of 20) for the purpose of this workflow. We filter to produce two subsets: a conservative set for species tree inference  (e.g. >= 4 species, mean copy \# <= 5), and one for which only gene family trees will be inferred (e.g. >= 4 species, mean copy \# <= 10).

We then distribute these two subsets of orthogroups into two, independent channels to ba anlyzed in parallel. Simultanously for both subsets, we infer multiple sequences alignments (using MAFFT), trim these alignments for gappy/uninformative regions (using ClipKit), and infer gene-family trees using IQ-TREE under a approximated mixture model of amino acid substitution (LG+C40+F+G) using the posterior mean site frequency model (PMSF).

Using the first, conservative subset of gene family trees, we infer a starting, unrooted species tree using Asteroid. This starting species tree is provided along with each gene family tree and corresponding multiple sequence alignment to SpeciesRax, which roots the species tree and improves the topology under a model of gene duplication, loss and transfer.

Using this improved species tree, we then use GeneRax for both subsets of gene families, reconciling them with the species tree and inferring rates of gene duplication, transfer and loss on a per-family and (TBD) per-species basis.

### The workflow thus proceeds as follows:
1. Proteomes are downloaded from S3
2. Each proteome is summarized using [`BUSCO`](https://busco.ezlab.org/) completeness using both Eukaryota and the taxon-specific level
3. Proteomes originating from UniProt (with UniProt protein accessions) are annotated using UniProt.ws
4. Proteomes are staged/reformated for analysis with OrthoFinder
5. Determine All-v-All sequence similarity using [`Diamond`](https://github.com/bbuchfink/diamond) BlastP ultra-sensitive
6. Cluster UniProt sequences into orthogroups/gene-families using [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)'s implementation of [`MCL`](http://micans.org/mcl/) clustering using a specified set of inflation scores
7. Summarization and quantification of orthogroup inference performance using a set of summary statistics, including the functional protein domain score using [`COGEQC`](https://almeidasilvaf.github.io/cogeqc/index.html)
8. Repeat step five (5: MCL clustering into orthogroups) for all species under the optimal inflation parameter
9. Summarize distribution of orthogroups across taxonomic groups and per-species copy number, filtering into a conservative subset for species tree inference, and one for gene-family tree inference.
10. Infer multiple sequence alignments for each focal gene family with [`MAFFT`](https://mafft.cbrc.jp/alignment/software/)
11. Trim uninformative/memory-consuming/gappy segments of alignments with [`ClipKit`](https://github.com/JLSteenwyk/ClipKIT)
12. Infer gene family trees using [`IQ-TREE`](http://www.iqtree.org/)
13. Infer a starting species tree using Asteroid [`Asteroid`](https://github.com/BenoitMorel/Asteroid)
14. Root the species tree, improving its topology under a model of gene duplication, transfer, and loss using [`SpeciesRax`](https://github.com/BenoitMorel/GeneRax)
15. Reconcile gene family trees with the species tree, inferring rates of gene duplication, transfer and loss using [`GeneRax`](https://github.com/BenoitMorel/GeneRax)
16. Infer phylogenetically hierarchical orthologs using [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)

#### Upon initiation, several things happen in parallel.
1. First, proteomes are summarized with BUSCO at (currently) two user-defined scales: one that is taxonomically shallow (as close to each species as possible), and one that is taxonomically broad (e.g. Eukaryotes).
   - The goal here is to quantify input proteome quality and completeness - comparison at the broad-scale (e.g. across eukaryotes) may facilitate more direct biological comparisons across species, but for understudied groups may provide misleadingly poor summaries of proteome completeness as compared to analysis at the shallower taxonomic scale.
   - Results of these analyses are independent of everything downstream and so each may be done in parallel of all other steps (i.e. will not prevent subsequent analyses from beginning).
2. Second, we prepare proteomes for the highly parallized distribution of all-v-all protein comparisons with Diamond BlastP using Nextflow and subsequent MCL clustering with OrthoFinder. This involves calling orthofinder, and formatting data in a manner that the software is built to deal with. This is largely a technicality.
3. Third, we begin running the all-v-all protein comparisons in parallel for both the MCL inflation parameter testing dataset and complete dataset, distributing jobs across AWS spot instances. Doing so in this way, we limit the extent to which the MCL inflation parameter test slows progress on analyzing the complete dataset, and takes advantage of the wealth of compute resources provided by AWS, all while saving up to 90% on costs by using spot instances.

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->


## Quick Start
TODO: ADD ORTHOGROUP FILTERS TO PARAMS.FILE AND PHYLORTHOLOGY.NF
1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)
2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Create two samplesheets (one for the complete dataset, one for MCL testing) following the required format:
NOTE: the "mode" column for busco can be incorporated into params file
```
species,file,taxonomy,shallow,broad,mode,uniprot,mcl_test
Entamoeba_histolytica,Entamoeba_histolytica-test-proteome.fasta,Amoebozoa,eukaryota_odb10,eukaryota_odb10,proteins,true,true
```
Description of Columns:  
`species`: species name to use. 
`file`: complete path to fasta file, whether local or remote (e.g. on S3 - provide S3 URI). 
`taxonomy`: higher-level taxonomy for each species (e.g. supergroup, class, family, genus). Utility of this parameter depends on the taxonomic scope of each dataset. Used in filtering orthogroups for phylogenetic inference.  
`shallow`: busco lineage dataset for shallow taxonomic scale analysis (e.g. below eukaryota). 
`broad`: busco lineage dataset for broad taxonomic scale analysis (e.g. eukaryota). 
`mode`: specification of busco analysis mode. 
`uniprot`: true/false specification indicating whether the proteome comes from UniProt (i.e. has UniProt protein accessions that we can use to annotate). 
`mcl_test`: true/false specification of whether this species is to be included in the MCL inflation parameter test-set. These species must have UniProt protein accessions (for COGEQC protein domain score).  

4. Create a parameter file that includes all necessary input, output, and parameter specifications - example below, followed by a description of parameters:
```
{
  "input": "/full/path/to/samplesheet.csv",
  "outdir": "results",
  "mcl_inflation": "1.0,2.0,3.0,4.0,5.0",
  "min_num_spp_per_og": 4,
  "min_num_grp_per_og": 2,
  "max_copy_num_spp_tree": 5,
  "max_copy_num_gene_trees": 10,
  "download_annots": "none",
  "tree_model": "LG+F+G4",
  "tree_model_pmsf": "LG+C40+F+G4"
}
```
`input`: Complete filepath to input samplesheet. May be locally stored, or remotely stored (e.g. on S3 - provide comple URI). 
`mcl_inflation`: Set of MCL inflation parameters to be tested when clustering proteins into orthogroups with OrthoFinder.  
`min_num_spp_per_og`: Minimum \# of species an orthogroup must contain for phylogenetic inference. 
`min_num_grp_per_og`: Minimum \# of 'higher' order taxonomic groups and orthogroup must contain for phylogenetic inference. 
`max_copy_num_spp_tree`: Maximum \# of per-species gene copy number an orthogroup may contain for species-tree inference. 
`max_copy_num_gene_trees`: Maximum \# of per-species gene copy number an orthogroup may contain for gene tree - species tree reconciliation with GeneRax.  
`download_annots`: Set of annotations to be downloaded. "none" corresponds to a minimal set. See description of parameters for expanded description of options.   
`tree_model`: Model of amino acid substition to be used for phylogenetic inference. If using a posterior mean site frequency model (see below), this model will be used to infer an initial guide-tree.  
`tree_model_pmsf`: OPTIONAL posterior mean site frequency model to be used for phylogenetic inference. If not specified (i.e. excluded from parameter file), only `tree_model` will be used.  


Alternatively, you can use the test dataset provided by Arcadia Science [here](https://github.com/Arcadia-Science/test-datasets/phylorthology).

5. Ensure that proteins are named following the following convention: `Species_genus:GeneID`  
```
# Example:  
>Entamoeba_histolytica:C4M6M9 tr|C4M6M9|C4M6M9_ENTHI NGG1-interacting factor 3, putative OS=Entamoeba histolytica HM-1:IMSS OX=294381 GN=EHI_161860 PE=3 SV=1
```
- Everything prior to the colon (:) is a constant (and unique to that species/lineage) identifier and can be whatever you would like, but must not include a colon. Everything that follows the colon is what must be a protein/gene identifier unique to that sequence. Additional sequence info may be included following a space. If the proteome comes from UniProt, this must be the UniProt protein accession. We use the string that follows the colon to extract the uniprot accession and annotate proteins.

6. Download the pipeline and test it on a minimal dataset with a single command run in the root of this repository:

   ```bash
   # If you're using your own test data
   nextflow run . -profile docker -params-file parameters.json
   ```

   OR

   ```bash
   # If you're using Arcadia Science's test dataset
   nextflow run . -profile docker -params-file https://github.com/Arcadia-Science/test-datasets/raw/main/phylorthology/nextflow_parameters.json
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

7. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run . -profile docker -params-file <PARAMS.JSON>
   ```
   
## Outputs
Below is a description of outputs, namely directory structure and their contents, produced from an end-to-end run of the workflow. This description will be updated as development progresses.  
1) `busco/`: Contains output of all BUSCO analyses. Primary directory includes:   
    - short summary results for each species, to each lineage dataset (shallow or broad) in text and json formats.  
    - batch summary results for each species, to each lineage dataset.  
    - directory for each species' busco analysis with more detailed output.  
2) `protein_annotations/`: protein annotations obtained per each species for which uniprot accessions (i.e. those corresponding to RefSeq protein accessions) are available.  
    - Which annotations are included here depends on how the "download_annots" parameter was specified.  
3) `diamond/`: contains the results of all pairwise comparisons of sequence similarity, between and within species using diamond BLASTP.  
    - TODO - would like you delete the TestBlast outputs at the conclusion of the workflow (or at least after their use/when running the complete mcl clustering.). - need to figure this out still.   
    - Can't be solved with the symlink, since it depends on whether we're MCL testing or not.  
4) `OrthoFinder/`: contains all results from orthofinder runs, including the MCL testing stage (`mcl_test_dataset`) and the full analysis using the best-performing inflation parameter (complete_dataset).  
    - `mcl_test_dataset/`: Contains one directory per tested inflation parameter. Directory structure follows OrthoFinder convention - orthogroup sequences are not retained to save space.  
        - TODO - delete persistent "OrthoFinder" directory. fixed with "publish as symlink"? Would be great if these links could just be ephemeral.  
    - `complete_dataset/`: Contains OrthoFinder output for the best-performing inflation parameter following standard convention (using OrthoFinder's -os flag).  
5) `orthogroup_summaries/`: Results from cogeqc.  
    - A table reporting all calculated summary statistics for the orthogroups inferred for each inflation parameter.  
    - A text file that lists the best-performing inflation parameter.  
    - A PDF plotting the summary statistics as a function of inflation parameter value.  
6) `filtered_orthogroups/`: final orthogroups used for gene-family tree and species tree inference that passed user-specifed filters.  
    - `species_tree_og_msas/`: trimmed multiple sequence alignments for gene families used in species tree inference.  
    - `gene_tree_og_msas/`: trimmed multiple sequence alignments for gene families used only in gene-family tree / species tree reconciliation.  
    - `all_ogs_counts.csv`: comma-separated csv listing, for all orthogroups (including those fro which msa/gene family trees are not inferred), the number of included species, total copy number, mean copy number, and number of higher-level taxonomic groups included.  
    - `(gene)speciestree_core_ogs_counts.csv`: the same as above, but for only the two respects subsets of gene families.  
7) `mafft/`: multiple sequence alignments for orthogroups passing filtering thresholds (e.g. minimum number of species, maximum copy number).  
8) `trimmed_msas/`: Multiple sequence alignments inferred using mafft, trimmed for phylogenetic inference with ClipKit.  
9) `iqtree/`: gene family trees (and corresponding log files) for all gene families for which (trimmed) multiple sequence alignments were estimates.  
10) `species_tree_prep/`: set of files used by asteroid and Gene/SpeciesRax to correspond gene-family protein sequences to species IDs, etc.  
11) `asteroid/`: starting, unrooted species tree inferred using all gene family trees with asteroid. Includes:  
    - asteroid.bestTree.newick: Single species tree with the greatest likelihood among set of inferred trees (for instance if multiple random starting trees are used).  
    - asteroid.allTrees.newick: All inferred species trees (=1 if default of 1 random starting tree is used).  
    - asteroid.scores.txt: likelihood scores for all inferred trees.  
12) `speciesrax/`: all results/outputs from SpeciesRax inferred using the subset of gene families that passed filters for involvement in rooted species tree inference. Full description of these outputs (including gene-family tree/species tree reconciliations) are described on the GeneRax github. Key output directories includes:  
    - `species_trees/`: contains inferred rooted species trees, species tree likihoods, and species trees with internal nodels labeled according to their support values.  
    - `reconciliations/`: gene-family tree - species tree reconciliations (species trees with gene family duplications, transfers, and losses mapped on) in several formats. May be plotted using reconciliation software like thirdkind. Additionally contains files for each gene family that enumerate duplication-transfer-loss event counts per species.  
    - `results/`: one directory per-gene-family containing reconciled gene trees and inferred rates of gene family duplication, transfer and loss.  
13) `generax/`: all results/outputs from GeneRax inferred using the subset of gene families that passed filters for involvement in rooted species tree inference. Full description of these outputs (including gene-family tree/species tree reconciliations) are described on the GeneRax github.  
    - Directory structure here is the same as for SpeciesRax, but there is no directory for species trees, as this step was conducted during the SpeciesRax module.  

## Credits

phylorthology was originally written by Arcadia Science.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  Arcadia-Science/phylorthology for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

## Parameter specification
TODO: This is a work in progress
We have set sensible parameter choices as default for each module, however several modules have parameters that are best-suited to user specification on a per-analysis basis. The section below describes, for each module, how these specifications are to be made as well as what some common specifications may be. Where necessary, refer to the documentation of each respective software for a more complete list of possible parameter choices.

### Module:
#### ANNOTATE_UNIPROT:
i) `download_annots`: Specified in the parameter file. Parameter may be specified as one of three things:  
    a) `all` - download all 16 possible sets of protein annotations from UniProt where possible.  
    b) `none` - download only the minimum necessary annotations that are used by cogeqc for orthogroup inference quality assessments.  
    c) A quoted, comma separated string of select numbers 1-16: example `"1,2,4,7,10"`. Numbers correspond to the index of annotations the user would like to download. See below for the correspondance and brief description of each annotation. For indices 4-16 (in particular) see https://www.uniprot.org/help/return_fields.  
        `1`) Minimal set of protein annotations/metadata required for COGEQC orthogroup inference. Include protein external IDs for InterPro, SupFam, ProSite, HOGENOM, OMA, and OrthoDB.  
        `2`) General protein metadata, including protein name, length, mass, information from mass spec experiments, host organisms (for viral proteins), which organelle (if relevant) encoding the protein, and any AA variants due to RNA editing.  
        `3`) Gene ontologies - biological process, cellular component, molecular function, ontology ID.  
        `4`) Function: Multiple annotations pertaining to the molecular function of the protein.  
        `5`) Interactions.  
        `6`) Protein-protein interactions (external database reference IDs).  
        `7`) Pathology & biotech.  
        `8`) Subcellular location.  
        `9`) Post-translation modification (PTM) / processsing.  
        `10`) PTM databases.  
        `11`) Protein family & domains.  
        `12`) Protein family/group databases.  
        `13`) Sequence databases.  
        `14`) 3D structure databases.  
        `15`) Enzyme and pathway databases.  
        `16`) Phylogenomic databases.  

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## TODOs before v0.0.1

- [ ] For each module in `nf-core-modified` document, git-sha, branch and repo of the module, along with what was adapted and why.
- [ ] Make sure the CONTRIBUTING doc is accurate and move it to the root of the repo
- [ ] Add a CODE_OF_CONDUCT doc
