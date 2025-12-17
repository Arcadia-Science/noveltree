# NovelTree: Highly parallelized phylogenomic inference

**Arcadia-Science/noveltree** is a phylogenomic pipeline designed to analyze proteomes from diverse organisms and inferring orthology, gene-family trees, and a species tree. The pipeline assumes that prior to analysis, input proteomes for each species have been sufficiently filtered such that no additional filtering of species or sequences is required. For a description of such a filtering procedure, see the following [GitHub repository](https://github.com/Arcadia-Science/2023-tsar-noveltree/tree/main/scripts/data-preprocessing).

![Workflow Figure](./Fig2-Workflow-part-one.png)
![Workflow Figure](./Fig4-Workflow-part-two.png)

*These figures illustrate the full workflow mode. Simplified and zoogle modes skip certain steps (e.g., BUSCO, per-family GeneRax) or add additional analyses (e.g., phylo-dist). See [Workflow Modes](#workflow-modes) for details.*

`NovelTree` is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

---

## Quick Start

**NOTE: Unfortunately, at this time NovelTree is not compatible with Apple silicon/ARM architectures (e.g. M1, M2 chips).**

**1.** Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`).

**2.** Install [`Docker`](https://docs.docker.com/engine/installation/).

**3.** Build the required Docker images (optional - only needed if you've modified the pipeline or are using a custom fork):

```bash
make docker-all
```

Or build individual images:
```bash
make docker-physicochemical-props
make docker-phylo-dist
```

**4.** Download the pipeline and our minimal test dataset with a single command run in the root of this repository:

```bash
nextflow run . -profile docker -params-file https://github.com/Arcadia-Science/test-datasets/raw/main/noveltree/tsar_downsamp_test_parameters.json
```

In cases where you need to specify resource usage limits to NovelTree (e.g. you are running it on a local desktop or laptop), you can specify the maximum available CPU and memory resources as follows:

```bash
nextflow run . -profile docker -params-file https://github.com/Arcadia-Science/test-datasets/raw/main/noveltree/tsar_downsamp_test_parameters.json  --max_cpus 12 --max_memory 16GB
```

Nextflow requires some memory resources to be allocated for overhead - consequently, we suggest reducing the specified `--max_memory` by ~2GB or more below the amount available to your particular compute environment.

**NOTE: The workflow supports both Docker and Singularity profiles.**

---

## Workflow Modes

NovelTree supports three workflow modes to accommodate different use cases and computational constraints:

| Feature | Full | Simplified | Zoogle |
|---------|:----:|:----------:|:------:|
| BUSCO quality assessment | ✓ | ✗ | ✗ |
| Default aligner | WITCH | FAMSA | FAMSA |
| Per-family GeneRax | ✓ | ✗ | ✗ |
| Per-species GeneRax | ✓ | ✓ | ✓ |
| GeneRax strategy | SPR | EVAL | EVAL |
| Phylogenetic profiles | ✓ | ✓ | ✓ |
| Physicochemical properties | ✗ | ✗ | ✓ |
| Time-calibrated species tree | ✗ | ✗ | ✓ |
| Phylo-dist analysis | ✗ | ✗ | ✓ |

### Full Mode (Default)

The complete pipeline with all optional analyses enabled. Best for comprehensive phylogenomic studies where accuracy is prioritized over speed.

```bash
nextflow run . -profile docker --input samplesheet.csv --outdir results
```

### Simplified Mode

A streamlined, high-throughput variant optimized for large datasets. Uses FAMSA (faster) instead of WITCH for alignment, skips BUSCO quality assessment, and runs only per-species GeneRax with the faster EVAL strategy.

```bash
nextflow run . -profile docker,simplified --input samplesheet.csv --outdir results
```

### Zoogle Mode

Inherits simplified mode settings and adds analyses for organism prioritization: physicochemical protein properties, time calibration of the species tree, and phylogenetically-corrected protein distance analysis. Requires a reference time-calibrated tree and specification of a reference species.

```bash
nextflow run . -profile docker,zoogle \
  --input samplesheet.csv \
  --outdir results \
  --reference_time_tree reference_timetree.newick \
  --ref_species Genus_species
```

---

## Running on AWS Batch

NovelTree includes a dedicated AWS Batch profile optimized for cloud-scale analyses:

```bash
nextflow run . \
  -profile awsbatch \
  --awsqueue <your-batch-queue> \
  --awsregion <your-aws-region> \
  -work-dir s3://<your-bucket>/work \
  --outdir s3://<your-bucket>/results \
  --input s3://<your-bucket>/samplesheet.csv
```

The `awsbatch` profile includes optimized executor settings (queue size of 1000 jobs) and automatic report overwriting for seamless pipeline resumption.

**Requirements:**
- AWS Batch compute environment and job queue configured
- Work directory (`-work-dir`) and output directory (`--outdir`) must be S3 paths
- Input samplesheet and proteome files accessible from S3
- Appropriate IAM permissions for Batch and S3 access

For detailed AWS Batch setup instructions, see the [usage documentation](docs/usage.md#running-on-aws-batch).

---

## Running with Singularity

NovelTree supports Singularity as an alternative to Docker, which is useful for HPC environments where Docker may not be available:

```bash
nextflow run . -profile singularity --input samplesheet.csv --outdir results
```

Docker images are automatically pulled and converted to Singularity format. Converted images are cached in `${outdir}/singularity_cache` to avoid repeated conversions on subsequent runs.

For detailed Singularity instructions, see the [Singularity documentation](docs/singularity.md).

---

## Building Docker Images

NovelTree uses custom Docker images for specific analysis modules. If you're using a custom fork or have modified the pipeline code, you'll need to rebuild these images.

### Quick Build

Build all required images:
```bash
make docker-all
```

### Individual Image Builds

Build specific images:
```bash
# Physicochemical properties calculation module
make docker-physicochemical-props

# Phylogenetic distance analysis module
make docker-phylo-dist
```

### Technical Details

All Docker images are built from the repository root with the build context set to ensure access to vendored code in `bin/raas/`. The images are built for `linux/amd64` platform for compatibility.

**Note:** Building R-based images (phylo-dist) may take 15-20 minutes due to package compilation.

### Vendored RAAS Code

The `bin/raas/` directory contains code vendored from the [raas-organism-prioritization](https://github.com/Arcadia-Science/raas-organism-prioritization) repository. See `bin/raas/README.md` for provenance details including source commit and modifications.

---

## Pipeline summary (Defaults)

At its core, `NovelTree` is a compilation of methods that facilitates user-customized phylogenomic inference from whole proteome amino acid sequence data. **_The method automates all steps of the process, from calculating reciprocal protein-sequence similarity to gene-family inference, multiple sequence alignment and trimming, gene-family and rooted species tree inference, to inference of gene-family evolutionary dynamics._**

Because `NovelTree` is built in [Nextflow](https://www.nextflow.io), the workflow distributes tasks in a highly parallel and asynchronous manner across available computational resources. The workflow supports multiple execution environments including local execution, [AWS Batch](#running-on-aws-batch) for cloud-scale analyses, and SLURM schedulers for HPC clusters ([see Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html)).

To account for the confounding effects of sequence length (and thus evolutionary) divergence on sequence similarity scores, `NovelTree` leverages [`OrthoFinder`](https://github.com/davidemms/OrthoFinder) to normalize these similarity scores prior to clustering into orthogroups/gene families with MCL clustering. Because this clustering is contingent upon the MCL inflation parameter, `NovelTree` automates the identification of the inflation parameter that returns the most biologically sensible set of orthogroups when a list of MCL inflation values is provided. If a single MCL inflation is provided by the user, the pipeline will use that as the best-performing inflation parameter. Based on our own [analyses](https://doi.org/10.57844/arcadia-z08x-v798), we would suggest using an inflation parameter of `2.5` if you elect to use a singular value.

**Thus, two rounds of protein clustering takes place when a list of MCL inflation parameters is provided. If a single MCL inflation parameter is provided by the user, the first step is skipped and the second step is run using the user-supplied inflation parameter as the best-performing one.:**
**1.** An initial round for inflation parameter testing on a (reduced) set of proteomes for which UniProt protein accessions are available, and
**2.** A second round on the complete dataset.

Once the first round of MCL clustering has completed, `NovelTree` summarizes orthogroups based on a number of metrics, choosing a best-performing inflation parameter for the analysis of the full dataset. This includes a functional protein annotation score calculated with [`COGEQC`](https://almeidasilvaf.github.io/cogeqc/index.html), which quantifies the ratio of InterPro domain "Homogeneity" within orthogroups to "Dispersal" of domains among orthogroups. This statistic is also calculated for OMA orthology database IDs.

With orthogroups/gene families inferred, `NovelTree` will summarize each gene family on the basis of their taxonomic and copy number distribution, quantifying the number of species/clades included in each, as well as the mean per-species copy number. These summaries facilitate 'filtering' for sufficiently conserved/computationally tractable gene families for downstream phylogenetic analysis. In other words, it may be best, depending on use-case, to avoid excessively small (e.g. < 4 species) or large gene families (e.g. > 50 species and mean copy # of 20 - this upper limit will depend on available computational resources) for the purpose of this workflow. We filter to produce two subsets: a conservative set for species tree inference (e.g. >= 4 species, mean copy \# <= 5), and one for which only gene family trees will be inferred (e.g. >= 4 species, mean copy \# <= 10).

For both subsets, `NovelTree` subsequently infers cleaned multiple sequences alignments (using [`WITCH`](https://github.com/c5shen/WITCH) by default, with [`MAFFT`](https://mafft.cbrc.jp/alignment/software/) or [`FAMSA`](https://github.com/refresh-bio/FAMSA) as alternatives) and gene-family trees using [`FastTree2`](http://www.microbesonline.org/fasttree/).

Using the first conservatively sized subset of gene family trees, `NovelTree` infers a starting, unrooted species tree using [`Asteroid`](https://github.com/BenoitMorel/Asteroid), a highly computationally efficient method. In parallel, a second species tree is inferred using [`SpeciesRax`](https://github.com/BenoitMorel/GeneRax/wiki/SpeciesRax), which roots the species tree reconciling the topology of the species tree with each gene family tree under a model of gene duplication, loss and transfer.

Using this improved species tree, `NovelTree` then uses [`GeneRax`](https://github.com/BenoitMorel/GeneRax) for both subsets of gene families, reconciling them with the species tree and inferring rates (and per-species event counts) of gene duplication, transfer and loss for each gene family and each species, using both the per-family, and per-species models.

With the rooted species tree inferred, `NovelTree` uses [`OrthoFinder`](https://github.com/davidemms/OrthoFinder) one final time to parse each orthogroup/gene family into phylogenetically hierarchical orthogroups.

### Phylogenetic Profiles

The GeneRax reconciliation outputs are curated into **phylogenetic profiles**—species × gene family matrices that summarize evolutionary events across the phylogeny. These matrices include:

- **Duplication counts**: Gene duplications per species-tree node per gene family
- **Loss counts**: Gene losses per species-tree node per gene family
- **Speciation counts**: Speciation events per species-tree node per gene family
- **HGT donor counts**: Horizontal gene transfer events where the species is the donor
- **HGT recipient counts**: Horizontal gene transfer events where the species is the recipient
- **HGT summed counts**: Combined donor and recipient transfer events

These profiles provide a comprehensive view of gene family evolutionary dynamics and can be used for downstream comparative analyses.

### Zoogle Mode Analyses

When running with the `zoogle` profile, NovelTree performs additional analyses designed for organism prioritization based on protein evolution:

1. **Physicochemical Properties**: Calculates amino acid composition and physicochemical properties (molecular weight, aromaticity, instability index, flexibility, hydropathy, isoelectric point, charge, and secondary structure fractions) for all proteins in each gene family.

2. **Time Calibration**: Calibrates the inferred species tree against a user-provided reference timetree using congruification, enabling evolutionary rate comparisons across lineages.

3. **Phylogenetically-Corrected Protein Distances**: For each gene family, computes multivariate distances between proteins based on their physicochemical properties, correcting for phylogenetic non-independence. Statistical tests identify proteins that are exceptionally (dis)similar to a reference species, which may be indicative of unusual evolutionary divergence, convergence, or conservatism.

---

## Usage
For a detailed description of basic- to advance-usage of the workflow, please see the [`usage.md`](docs/usage.md) file.

## Outputs
For a detailed description of workflow outputs, please see the [`outputs.md`](docs/outputs.md) file.

---

## Credits

NovelTree was originally written by Arcadia Science.

## Feedback, contributions, and reuse

We try to be as open as possible with our work and make all of our code both available and usable.
We love receiving feedback at any level, through comments on our pubs or Twitter and issues or pull requests here on GitHub.
In turn, we routinely provide public feedback on other people’s work by [commenting on preprints](https://sciety.org/lists/f8459240-f79c-4bb2-bb55-b43eae25e4f6), filing issues on repositories when we encounter bugs, and contributing to open-source projects through pull requests and code review.

Anyone is welcome to contribute to our code.
When we publish new versions of pubs, we include a link to the "Contributions" page for the relevant GitHub repo in the Acknowledgements/Contributors section.
If someone’s contribution has a substantial impact on our scientific direction, the biological result of a project, or the functionality of our code, the pub’s point person may add that person as a formal contributor to the pub with "Critical Feedback" specified as their role.

Our policy is that external contributors cannot be byline-level authors on pubs, simply because we need to ensure that our byline authors are accountable for the quality and integrity of our work, and we must be able to enforce quick turnaround times for internal pub review.
We apply this same policy to feedback on the text and other non-code content in pubs.

If you make a substantial contribution, you are welcome to publish it or use it in your own work (in accordance with the license — our pubs are CC BY 4.0 and our code is openly licensed).
We encourage anyone to build upon our efforts.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
If you use  Arcadia-Science/noveltree for your analysis, please cite it using the following doi: [10.57844/arcadia-z08x-v798](https://doi.org/10.57844/arcadia-z08x-v798)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

---

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
