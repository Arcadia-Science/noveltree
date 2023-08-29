# NovelTree: Highly parallelized phylogenomic inference  

**Arcadia-Science/noveltree** is a phylogenomic pipeline designed to analyze proteomes from diverse organisms and inferring orthology, gene-family trees, and a species tree. The pipeline assumes that prior to analysis, input proteomes for each species have been sufficiently filtered such that no additional filtering of species or sequences is required. For a description of such a filtering procedure, see the following [GitHub repository](https://github.com/Arcadia-Science/2023-tsar-noveltree/tree/main/scripts/data-preprocessing).

![Workflow Figure](https://github.com/Arcadia-Science/noveltree/blob/main/NovelTree.png)  

`NovelTree` is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.  

---

## Quick Start

**1.** Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`).  

**2.** Install [`Docker`](https://docs.docker.com/engine/installation/).  

**3.** Download the pipeline and our minimal test dataset with a single command run in the root of this repository:  

```bash
nextflow run . -profile docker -params-file https://github.com/Arcadia-Science/test-datasets/raw/main/noveltree/tsar_downsamp_test_parameters.json
```

**NOTE: Currently the workflow only works using the docker profile.**  

---

## Pipeline summary (Defaults)

At its core, `NovelTree` is a compilation of methods that facilitates user-customized phylogenomic inference from whole proteome amino acid sequence data. **_The method automates all steps of the process, from calculating reciprocal protein-sequence similarity to gene-family inference, multiple sequence alignment and trimming, gene-family and rooted species tree inference, to inference of gene-family evolutionary dynamics._**  

Because `NovelTree` is built in [Nextflow](https://www.nextflow.io), the workflow distributes tasks in a highly parallel and asynchronous manner across available computational resources. The workflow is currently optimized for a single computational environment but is continually being developed for deployment across AWS spot-instances with Nextflow Tower, and may also be configured to run in a highly parallel manner on SLURM schedulers ([see here for documentation](https://www.nextflow.io/docs/latest/executor.html)).  

To account for the confounding effects of sequence length (and thus evolutionary) divergence on sequence similarity scores, `NovelTree` leverages [`OrthoFinder`](https://github.com/davidemms/OrthoFinder) to normalize these similarity scores prior to clustering into orthogroups/gene families with MCL clustering. Because this clustering is contingent upon the MCL inflation parameter, `NovelTree` automates the identification of the inflation parameter that returns the most biologically sensible set of orthogroups.  

**Thus, two rounds of protein clustering takes place:**  
**1.** An initial round for inflation parameter testing on a (reduced) set of proteomes for which UniProt protein accessions are available, and  
**2.** A second round on the complete dataset.  

Once the first round of MCL clustering has completed, `NovelTree` summarizes orthogroups based on a number of metrics, choosing a best-performing inflation parameter for the analysis of the full dataset. This includes a functional protein annotation score calculated with [`COGEQC`](https://almeidasilvaf.github.io/cogeqc/index.html), which quantifies the ratio of InterPro domain "Homogeneity" of domains within orthogroups to "Dispersal" of domains among orthogroups. This statistic is also calculated for OMA orthology database IDs.  

With orthogroups/gene families inferred, `NovelTree` will summarize each gene family on the basis of their taxonomic and copy number distribution, quantifying the number of species/clades included in each, as well as the mean per-species copy number. These summaries facilitate 'filtering' for sufficiently conserved/computationally tractable gene families for downstream phylogenetic analysis. In other words, it may be best, depending on use-case, to avoid excessively small (e.g. < 4 species) or large gene families (e.g. > 50 species and mean copy # of 20 - this upper limit will depend on available computational resources) for the purpose of this workflow. We filter to produce two subsets: a conservative set for species tree inference (e.g. >= 4 species, mean copy \# <= 5), and one for which only gene family trees will be inferred (e.g. >= 4 species, mean copy \# <= 10).  

For both subsets, `NovelTree` subsequently infers cleaned multiple sequences alignments (using [`WITCH`](https://github.com/c5shen/WITCH)) and gene-family trees using [`FastTree2`](http://www.microbesonline.org/fasttree/).  

Using the first conservatively sized subset of gene family trees, `NovelTree` infers a starting, unrooted species tree using [`Asteroid`](https://github.com/BenoitMorel/Asteroid), a highly computationally efficient method. In parallel, a second species tree is inferred using [`SpeciesRax`](https://github.com/BenoitMorel/GeneRax/wiki/SpeciesRax), which roots the species tree reconciling the topology of the species tree with each gene family tree under a model of gene duplication, loss and transfer.  

Using this improved species tree, `NovelTree` then uses [`GeneRax`](https://github.com/BenoitMorel/GeneRax) for both subsets of gene families, reconciling them with the species tree and inferring rates (and per-species event counts) of gene duplication, transfer and loss for each gene family and each species, using both the per-family, and per-species models.  

With the rooted species tree inferred, `NovelTree` uses [`OrthoFinder`](https://github.com/davidemms/OrthoFinder) one final time to parse each orthogroup/gene family into phylogenetically hierarchical orthogroups.  

---

## Usage
For a detailed description of basic- to advance-usage of the workflow, please see the [`usage.md`](docs/usage.md) file. 

## Outputs
For a detailed description of basic- to advance-usage of the workflow, please see the [`outputs.md`](docs/outputs.md) file. 

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
<!-- If you use  Arcadia-Science/noveltree for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

---

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## TODOs before v0.0.1

- [ ] For each module in `nf-core-modified` document, git-sha, branch and repo of the module, along with what was adapted and why.
- [ ] Make sure the CONTRIBUTING doc is accurate and move it to the root of the repo
- [ ] Add a CODE_OF_CONDUCT doc
