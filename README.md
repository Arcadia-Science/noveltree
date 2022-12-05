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
species,file,taxonomy,shallow,broad,mode,uniprot,mcl_test
Entamoeba_histolytica,Entamoeba_histolytica-test-proteome.fasta,Amoebozoa,eukaryota_odb10,eukaryota_odb10,proteins,true,true
- "uniprot" column is a true/false specification indicating whether the proteome comes from UniProt (i.e. has UniProt protein accessions that we can use to annotate)
- mcl_test is a specification of whether this species is to be included in the MCL inflation parameter test-set. These species must have UniProt protein accessions (for COGEQC protein domain score). 

4. Create a parameter file that includes all necessary input, output, and parameter specifications: example below:
{
  "input": "/path/to/complete_dataset_samplesheet.csv",
  "mcl_test_input": "/path/to/mcl_testing_dataset_samplesheet.csv",
  "s3_dir": "s3://path/to/aws-S3-directory/containing/fastas",
  "outdir": "name_of_results_dir",
  "mcl_inflation": [ "1.0", "1.1", "1.2" ]
}

5. Ensure that proteins are named following the following convention. 
- Species_genus:GeneID
- Example: >Entamoeba_histolytica:C4M6M9 tr|C4M6M9|C4M6M9_ENTHI NGG1-interacting factor 3, putative OS=Entamoeba histolytica HM-1:IMSS OX=294381 GN=EHI_161860 PE=3 SV=1
- Everything prior to the colon (:) is a constant (and unique to that species/lineage) identifier and can be whatever you would like, but must not include a colon. 
- Everything that follows the colon is what must be a protein/gene identifier unique to that sequence. 
- Additional sequence info may be included following a space. 
- If the proteome comes from UniProt, this must be the UniProt protein accession.
- We use the string that follows the colon to extract the uniprot accession and annotate proteins. 

6. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run Arcadia-Science/phylorthology -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run Arcadia-Science/phylorthology --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Credits

Arcadia-Science/phylorthology was originally written by Arcadia Science.

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

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
