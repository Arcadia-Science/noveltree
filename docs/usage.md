# Run the workflow on your own custom dataset

NovelTree includes an optional built-in preprocessing step (`--preprocess`) that handles TransDecoder for transcriptomes, isoform filtering, minimum protein length filtering, and redundancy removal. When enabled, per-species flags in the samplesheet control which steps run for each species.

Alternatively, the pipeline assumes that input proteomes have been sufficiently filtered externally. For a description of such a filtering procedure, see the following [GitHub repository](https://github.com/Arcadia-Science/2023-tsar-noveltree/tree/main/scripts/data-preprocessing).

## Verify your setup

```bash
nextflow run . -profile docker,test --outdir test_results
```

## Preparation

**1.** Prepare a samplesheet following the required format:

```
species,input_data,input_type
Homo-sapiens,UP000005640,proteins
```

The first 3 columns are required. All remaining columns are optional — include any subset in any order.

### Column Reference

| Column | Required | Default | Allowed Values | Description |
|--------|----------|---------|----------------|-------------|
| `species` | yes | — | `Genus-species` | Species name. Spaces and underscores are auto-converted to hyphens |
| `input_data` | yes | — | See [Data Sources](#data-sources) | Protein FASTA, nucleotide FASTA, URL, NCBI accession, or UniProt proteome ID |
| `input_type` | yes | — | `proteins` / `transcriptome` | `proteins` for amino acid input; `transcriptome` for nucleotide input |
| `has_uniprot_ids` | no | `no` | `yes` / `no` | Whether FASTA headers contain UniProt accessions (needed for MCL testing and reference proteome validation) |
| `include_in_mcl_test` | no | `no` | `yes` / `no` | Include in MCL inflation parameter testing. Requires `has_uniprot_ids=yes` |
| `transdecoder` | no | `no` | `yes` / `no` | Run TransDecoder ORF prediction. Only valid with `input_type=transcriptome`. Auto-disabled if `input_type=proteins` |
| `filter_isoforms` | no | `no` | `yes` / `no` | Filter to longest isoform per gene. Auto-set to `yes` for NCBI accessions |
| `reference_proteome` | no | `no` | `yes` / `no` | Curated reference proteome — skips redundancy removal. Requires `has_uniprot_ids=yes` |
| `busco_shallow` | no | `NA` | BUSCO lineage or `NA` | e.g. `primates_odb10`. Only used with `--busco true` |
| `busco_broad` | no | `NA` | BUSCO lineage or `NA` | e.g. `eukaryota_odb10`. Only used with `--busco true` |

### Data Sources

How to fill the `input_data` column:

**a) UniProt proteome ID** — e.g. `UP000005640`

Preferred for model organisms with curated reference proteomes. The pipeline queries the UniProt Proteomes API, resolves the taxonomy, and downloads the one-protein-per-gene FASTA from UniProt's FTP server.

- **When to use:** The organism has a reference proteome on [UniProt Proteomes](https://www.uniprot.org/proteomes/)
- **Format:** `UP` followed by 9+ digits
- **Recommended flags:** `input_type=proteins`, `has_uniprot_ids=yes`, `reference_proteome=yes`
- **Finding IDs:** Search UniProt Proteomes for your organism. Look for entries marked "Reference proteome"

**b) NCBI genome accession** — e.g. `GCA_030068145.1` or `GCF_000240725.1`

For organisms with annotated genomes on NCBI. The pipeline uses the NCBI `datasets` CLI to download and extract the `protein.faa` file. Isoform filtering is auto-enabled since NCBI proteomes typically include multiple isoforms per gene.

- **When to use:** The organism has an annotated genome on NCBI with protein predictions
- **Format:** `GCA_` or `GCF_` followed by digits, optionally with version (e.g. `.1`)
- **Recommended flags:** `input_type=proteins` (isoform filtering is auto-set to `yes`)
- **Verifying annotations:** Check the genome page on NCBI — the "Protein-coding genes" count should be > 0

**c) Direct URL to protein FASTA** — full URL to `.fasta.gz`, `.fa.gz`, etc.

For explicit control over which file to download. Works with any HTTP/HTTPS/FTP/S3 URL pointing to a protein FASTA.

- **When to use:** You have a specific URL to a protein FASTA file
- **Recommended flags:** `input_type=proteins`, `filter_isoforms=yes` if the source includes multiple isoforms

**d) Direct URL to nucleotide transcriptome** — full URL to `.fsa_nt.gz` or `.fna.gz`

For organisms with only transcriptome data (e.g., Transcriptome Shotgun Assembly). The pipeline downloads the file and runs TransDecoder to predict ORFs.

- **When to use:** The organism has a transcriptome assembly but no annotated genome
- **Format:** URL ending in `.fsa_nt.gz` or `.fna.gz`
- **Recommended flags:** `input_type=transcriptome`, `transdecoder=yes`, `filter_isoforms=yes`
- **Warning:** TransDecoder is designed for transcripts, NOT genome assemblies. Introns in genomic sequences will break ORF prediction

**e) Local file path** — absolute or relative path to a FASTA file

For files already on disk. Set flags based on the content type (protein vs. nucleotide).

- **Accepted extensions:** `.fasta`, `.fa`, `.fasta.gz`, `.fa.gz` (protein); additionally `.fna`, `.fna.gz` (nucleotide, with `transdecoder=yes`)

### Choosing a Data Source

```
Does the organism have a UniProt reference proteome?
  YES → Use UniProt ID, input_type=proteins, reference_proteome=yes
  NO  ↓
Does it have an annotated genome on NCBI (protein.faa available)?
  YES → Use NCBI accession, input_type=proteins
  NO  ↓
Does it have a Transcriptome Shotgun Assembly (TSA)?
  YES → Use URL to .fsa_nt.gz, input_type=transcriptome, transdecoder=yes
  NO  ↓
Do you have a local FASTA (protein or nucleotide)?
  → Use local path, set input_type accordingly
```

### Example Samplesheet

```csv
species,input_data,input_type,has_uniprot_ids,include_in_mcl_test,filter_isoforms,reference_proteome,transdecoder
Homo-sapiens,UP000005640,proteins,yes,yes,no,yes,no
Danio-rerio,GCF_000002035.6,proteins,no,no,yes,no,no
Nannochloropsis-sp,GCF_000240725.1,proteins,no,no,yes,no,no
Schizosaccharomyces-pombe,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/945/GCF_000002945.1_ASM294v2/GCF_000002945.1_ASM294v2_rna_from_genomic.fna.gz,transcriptome,no,no,yes,no,yes
Saccharomyces-cerevisiae,UP000002311,proteins,yes,no,yes,yes,no
Neurospora-crassa,euk_test_data/proteomes/Neurospora-crassa.fasta,proteins,yes,no,no,no,no
```

| Species | Source Type | Notes |
|---------|-------------|-------|
| Homo sapiens | UniProt ID | Reference proteome, MCL testing enabled |
| Danio rerio | NCBI accession | Annotated genome, isoform filtering |
| Nannochloropsis sp. | NCBI accession | Annotated genome, isoform filtering |
| S. pombe | URL (nucleotide) | Transcriptome, TransDecoder + isoform filtering |
| S. cerevisiae | UniProt ID | Reference proteome |
| N. crassa | Local file | Pre-existing protein FASTA |

### Preprocessing Behavior (when `--preprocess true`)

| `reference_proteome` | `transdecoder` | `filter_isoforms` | What Happens |
|---------------------|---------------|-------------------|--------------|
| `yes` | `no` | any | Quality cleanup only (no CD-HIT) |
| `no` | `yes` | `yes` | TransDecoder → isoform filter → CD-HIT 97% |
| `no` | `no` | `yes` | Isoform filter → CD-HIT 100% (exact dedup) |
| `no` | `no` | `no` | Quality cleanup → CD-HIT 100% (exact dedup) |

Quality cleanup always runs: removes stop codons, replaces rare amino acids (U→C, J/B/Z→X), and filters sequences shorter than `--min_protein_length` (default: 50 aa).

### FASTA Headers

When using NCBI/UniProt sources or `--preprocess`, the pipeline's `RENAME_FASTAS` module auto-standardizes headers to `Species-name_ProteinID`. Manual header formatting is only needed for local files used without `--preprocess`.

For local files without `--preprocess`, proteins must follow this convention: `Species-name:ProteinID`. If `has_uniprot_ids=yes`, the protein ID must be a UniProt accession.

**2.** Create a parameter file (see [Parameters](#parameters) for all options):

```json
{
  "input": "/full/path/to/samplesheet.csv",
  "outdir": "results",
  "mcl_inflation": "1.5",
  "aligner": "adaptive",
  "msa_trimmer": "clipkit",
  "tree_method": "iqtree",
  "ref_species": "Genus-species"
}
```

**3.** Start running your own analysis!

```bash
nextflow run . -profile docker -params-file <PARAMS.JSON>
```

Alternatively, you can use the test dataset provided by Arcadia Science [here](https://github.com/Arcadia-Science/test-datasets/noveltree).

---

## Parameters

### Core

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input` | — | Complete filepath to input samplesheet (local path, S3 URI, or URL) |
| `outdir` | `results` | Output directory path |
| `ref_species` | `Homo-sapiens` | Reference species for distance comparisons. Set to `'none'` for centroid-only analysis. Must match a samplesheet species name (format: `Genus-species`). Underscores and spaces are auto-converted to hyphens |
| `ncbi_email` | — | Email for NCBI Entrez queries. Required for zoogle mode when `reference_time_tree` is not provided |
| `reference_time_tree` | — | Path to a reference time-calibrated tree (Newick). If not provided in zoogle mode, a chronogram is auto-built from TimeTree.org |
| `test_run` | `false` | Restrict analysis to gene families containing all species for fast smoke testing |
| `test_run_mcl` | `false` | Enable MCL inflation parameter testing on species with UniProt accessions. When false, uses `--mcl_inflation` directly |
| `preprocess` | `false` | Enable built-in proteome preprocessing (TransDecoder, isoform filtering, min length, redundancy removal) |
| `simplified` | `true` | Enable simplified mode (skip BUSCO, per-species GeneRax EVAL only) |
| `busco` | `false` | Enable BUSCO quality assessment (full mode only) |
| `generax_per_family` | `false` | Enable per-family GeneRax analysis (full mode only) |

### Orthology & Filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mcl_inflation` | `"1.5"` | Quoted, comma-separated MCL inflation parameters for orthogroup clustering. A single value skips testing. Recommended single value: `2.5` ([rationale](https://doi.org/10.57844/arcadia-z08x-v798)) |
| `min_num_spp_per_og` | `4` | Minimum number of species a gene family must contain for phylogenetic inference |
| `min_prop_spp_for_spptree` | `0.50` | Minimum proportion of species for inclusion in species tree inference |
| `max_copy_num_spp_tree` | `10` | Maximum per-species gene copy number for species tree inference |
| `speciesrax_bundle_size` | `128` | Numeric orthogroup IDs per uncompressed SpeciesRax tree/mapping staging shard |
| `speciesrax_min_species_occupancy` | `0.50` | Temporary SpeciesRax-local minimum fraction of represented species |
| `speciesrax_max_mean_copies` | `4.0` | Temporary SpeciesRax-local maximum mean copies among represented species |
| `speciesrax_max_copies_per_species` | `8` | Temporary SpeciesRax-local maximum copies in any one species |
| `speciesrax_max_total_leaves_factor` | `4.0` | Temporary SpeciesRax-local maximum total leaves as a multiple of the species count |
| `min_protein_length` | `50` | Minimum amino acid sequence length during preprocessing (only when `--preprocess` enabled) |

### Alignment

| Parameter | Default | Description |
|-----------|---------|-------------|
| `aligner` | `"adaptive"` | Alignment method. Options: `adaptive` (MAFFT → WITCH → FAMSA by family size), `witch`, `mafft`, `famsa` |
| `align_tier1_max` | `200` | Max sequences for tier-1 MAFFT alignment; families above this use WITCH |
| `align_tier2_max` | `1000` | Max sequences for tier-2 WITCH alignment; families above this use FAMSA |
| `mafft_mode` | `"einsi"` | MAFFT algorithm for tier-1. Options: `einsi` (E-INS-i, conserved motifs) or `linsi` (L-INS-i, globally alignable) |
| `msa_trimmer` | `"clipkit"` | Trimming method. Options: `clipkit`, `cialign`, `none` |
| `min_ungapped_length` | `50` | Minimum ungapped length of cleaned/trimmed alignments |

### Tree Inference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `tree_method` | `"iqtree"` | Tree inference method. Options: `iqtree`, `fasttree`. IQ-TREE failures automatically fall back to FastTree |
| `tree_iqtree_max` | `1000` | Maximum sequences for IQ-TREE; larger families go directly to FastTree |
| `tree_model` | `"LG+F+G4"` | Amino acid substitution model for phylogenetic inference |
| `outgroups` | — | Comma-separated species IDs for manual rooting of the Asteroid species tree. If set, SpeciesRax estimates branch lengths on this tree |
| `iqtree_fasttree_fallback` | `true` | Automatically fall back to FastTree for IQ-TREE failures |

### Time Calibration (Zoogle Mode)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `time_calibration_method` | `"treePL"` | Species tree calibration method. Options: `treePL`, `PATHd8`. Gene family trees always use treePL with fixed smoothing |
| `age_bracket` | `0.20` | Proportional bracket applied to calibration ages for both species tree and gene family tree dating (e.g., 0.20 = ±20%) |

---

## Workflow Modes

See the [README Workflow Modes section](../README.md#workflow-modes) for full descriptions and command examples for Full, Simplified, and Zoogle modes.

### Additional examples

**Centroid-only mode** (no reference species):

```bash
nextflow run . -profile docker,zoogle \
  -params-file params.json \
  --ncbi_email user@example.com \
  --ref_species none
```

When `--ref_species none` is set, only centroid-based analyses are produced. Every gene family meeting the minimum size requirements (≥4 proteins, ≥2 species with ≥2 proteins each) will be analyzed.

**Combining profiles** with container engines:

```bash
# Simplified mode with Singularity
nextflow run . -profile singularity,simplified -params-file params.json

# Zoogle mode on AWS Batch (auto-build chronogram)
nextflow run . -profile awsbatch,zoogle \
  --awsqueue my-queue \
  --awsregion us-east-1 \
  -work-dir s3://bucket/work \
  --outdir s3://bucket/results \
  --ncbi_email user@example.com \
  --ref_species Genus-species
```

> **Supported container engines:** NovelTree is designed and tested to run with **Docker** or **Singularity/Apptainer**. While `conda`/`mamba` profiles exist (from the nf-core template), NovelTree does not ship curated Conda environments for its custom modules, so these profiles are not supported.

---

## Running on AWS Batch / Singularity

See the README for [AWS Batch](../README.md#running-on-aws-batch) and [Singularity](../README.md#running-with-singularity) setup instructions.

### Nextflow Tower (Publication Example)

When applying NovelTree to the dataset used in [the associated pub](https://doi.org/10.57844/arcadia-z08x-v798), we launched the workflow via Nextflow Tower to run on AWS Batch with:

1. `max_cpus = 5000`: Maximum CPUs (as spot instances) across all concurrent processes.
2. `max_memory = 30000.GB`: Maximum memory across all concurrent processes.
3. `max_time = 2400.h`: Maximum time across all concurrent processes.
4. 32 CPUs allocated to the head node for efficient job monitoring and submission.

---

## Pipeline steps

1. `INPUT_CHECK`: Proteomes are staged locally (including downloaded from S3 or other cloud storage if necessary)
2. When a list of mcl inflation values is provided, the pipeline performs these additional steps to select the best-performing MCL inflation parameter on a reduced set of proteomes for which UniProt protein accessions are available:
   1. `PROTEIN_ANNOTATION`: Proteomes for which sequence names include [`UniProt`](https://www.uniprot.org/) protein accessions are annotated with InterPro domains via the [UniProt ID Mapping API](https://www.uniprot.org/help/id_mapping)
   2. `ORTHOFINDER_PREP`: Proteomes are staged/reformated for analysis with [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)
   3. `DIAMOND_BLASTP`: Determine all-v-all (within and among species) protein sequence similarity using [`Diamond`](https://github.com/bbuchfink/diamond) BlastP ultra-sensitive
   4. `ORTHOFINDER_MCL`: Cluster [`UniProt`](https://www.uniprot.org/) sequences into orthogroups/gene-families using [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)'s implementation of [`MCL`](http://micans.org/mcl/) clustering using a specified set of inflation scores
   5. `COGEQC`: Summarization and quantification of gene family inference performance using a set of summary statistics, including the functional annotation score using [`COGEQC`](https://almeidasilvaf.github.io/cogeqc/index.html) applied to [`InterPro`](https://ebi.ac.uk/interpro/) domain annotations.
   6. `SELECT_INFLATION`: Based on the above summaries, select the (mean) inflation parameter that performs best (e.g. orthogroups are most homogenous in protein domain annotations, penalizing against dispersal of annotations across orthogroups), accounting for diminishing returns with increasing or decreasing parameter values.
3. `BUSCO` _(full mode only)_: Each proteome is summarized using [`BUSCO`](https://busco.ezlab.org/) completeness at both user-specified shallow (e.g. Eukaryota) and taxon-specific scales
4. `ORTHOFINDER_PREP`: All proteomes are staged/reformated for analysis with [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)
5. `DIAMOND_BLASTP`: Determine all-v-all (within and among species) protein sequence similarity using [`Diamond`](https://github.com/bbuchfink/diamond) BlastP ultra-sensitive.
6. `ORTHOFINDER_MCL`: Cluster [`UniProt`](https://www.uniprot.org/) sequences into orthogroups/gene-families using [`OrthoFinder`](https://github.com/davidemms/OrthoFinder)'s implementation of [`MCL`](http://micans.org/mcl/) clustering using a specified set of inflation scores
7. `ORTHOFINDER_MCL` also summarizes the distribution of orthogroups across species and per-species copy number, filtering into a conservative subset for species tree inference and one for gene-family tree inference.



8. `ALIGN_SEQS`: Infer multiple sequence alignments for each focal gene family using the adaptive three-tier strategy ([`MAFFT`](https://mafft.cbrc.jp/alignment/software/) for ≤200 seqs, [`WITCH`](https://github.com/c5shen/WITCH) for ≤3000, [`FAMSA`](https://github.com/refresh-bio/FAMSA) for larger), or a single aligner if specified
9. `TRIM_SEQS` _(optional)_: Trim uninformative/memory-consuming/gappy segments of alignments with either [`CIAlign`](https://github.com/KatyBrown/CIAlign) or [`ClipKit`](https://jlsteenwyk.com/ClipKIT/)
10. `INFER_TREES`: Infer gene family trees using either [`IQ-TREE`](http://www.iqtree.org/) (default) or [`FastTree2`](http://www.microbesonline.org/fasttree/)
11. `BUILD_REFERENCE_CHRONOGRAM`: When no explicit outgroup is provided, load a user-supplied rooted chronogram or build one from TimeTree.org. Zoogle mode reuses this exact tree for calibration.
12. `SPECIESRAX`: Infer the MiniNJ topology under MPI, transfer the encoded root split of the trusted chronogram (or explicit outgroups) onto that exact topology, then estimate branch lengths and quartet support under a duplication/loss model with [`SpeciesRax`](https://github.com/BenoitMorel/GeneRax/wiki/SpeciesRax). The final pass uses `--si-strategy SKIP`, so SpeciesRax preserves the transferred root. Core trees and mappings are validated and staged in bounded uncompressed shards, and alignments are omitted because gene-tree optimization is disabled.
13. `GENERAX_PER_FAMILY` _(full mode only)_: Reconcile gene family trees with the species tree, inferring rates of gene duplication and loss using [`GeneRax`](https://github.com/BenoitMorel/GeneRax) under the per-family model (rates are constant across all species/branches)
14. `GENERAX_PER_SPECIES`: Reconcile gene family trees with the species tree, inferring rates of gene duplication and loss using [`GeneRax`](https://github.com/BenoitMorel/GeneRax) under the per-species model (each species/branch has own rates). Uses SPR strategy in full mode, EVAL strategy in simplified/zoogle modes.
15. `PARSE_PHYLOHOGS`: Parse ortholog/paralog relationships and HOG membership from GeneRax reconciliation output
16. `PHYLO_PROFILES`: Generate phylogenetic profiles from each family's small GeneRax event and coverage files, summarizing gene duplication, loss, and speciation events across species and gene families
17. `TIME_CALIBRATE_SPECIES_TREE` _(zoogle mode only)_: Validate that the SpeciesRax root still matches the shared reference chronogram, then time-calibrate it using treePL penalized likelihood
18. `DATE_GENE_FAMILY_TREES` _(zoogle mode only)_: Time-calibrate gene family trees using speciation node ages from the dated species tree. Only speciation nodes from GeneRax reconciliation are used as calibration points.
19. `PROTEIN_PROPERTIES` _(zoogle mode only)_: Calculate amino acid composition and physicochemical properties for all gene families
20. `ZOOGLE` _(zoogle mode only)_: Calculate phylogenetically-corrected protein distances using Mahalanobis distances and permutation tests



---

# Advanced Usage

The sections below describe advanced usage of `NovelTree`, including per-module parameter specifications and outputs.

## Per-module parameter specification

We have set sensible parameter choices as default for each module, however several modules have parameters that are best-suited to user specification on a per-analysis basis. This section describes, for each module, fixed parameter names, or what default parameters specifications may be. Where necessary, refer to the documentation of each respective software for a more complete list of possible parameter choices.

Certain modules have parameters/flags that are specified in [`conf/modules.config`](../conf/modules.config); these are indicated as necessary. Custom specifications may be made following the same convention (example below, [documented here](https://nf-co.re/developers/modules#general)):

```
process {
    withName: 'MAFFT' {
        ext.args = [
            '--localpair',
            '--maxiterate 1000',
            '--anysymbol'
        ].join(' ')
    }
}
```

### Module: select each to follow links to corresponding module file.

#### 1. [`BUSCO`](../modules/nf-core-modified/busco.nf):

- `config_file`: Optional config file used used by BUSCO.
- `busco_lineages_path`: Optional path to locally stored BUSCO lineage datasets
- [BUSCO documentation](https://busco.ezlab.org/busco_userguide.html)

#### 2. [`ANNOTATE_UNIPROT`](../modules/local/annotate_uniprot.nf):

- Downloads InterPro domain annotations required for COGEQC gene family quality assessment via the UniProt ID Mapping API.

#### 3. [`DIAMOND_BLASTP`](../modules/nf-core-modified/diamond_blastp.nf):

- `--ultra-sensitive`: Specified in [`conf/modules.config`](../conf/modules.config). By default, sequence similarity is assessed using the most sensitive (albeit slowest) method.
- [Diamond documentation](https://github.com/bbuchfink/diamond/wiki)

#### 4. [`ORTHOFINDER_MCL`](../modules/local/orthofinder_mcl.nf):

- `mcl_inflation`: Comma-separated list of inflation parameter values to be used in testing. MCL inflation testing is optional. Set `--test_run_mcl true` to enable. When disabled (default), the pipeline uses the `--mcl_inflation` value directly.
- OrthoFinder's native membership and gene-count tables remain unchanged.
  Normal Nextflow task caching and stable published outputs provide `-resume`;
  NovelTree does not create a separate OrthoFinder checkpoint archive.
- [OrthoFinder2 documentation](https://github.com/davidemms/OrthoFinder)

#### 5. `ALIGN_SEQS`

##### [`MAFFT`](../modules/nf-core-modified/mafft.nf):

- Parameters specified in [`conf/modules.config`](../conf/modules.config). See MAFFT documentation for detailed description of options.
- `--localpair --maxiterate 1000 --anysymbol`: Runs MAFFT L-INS-i. Iterative refinement method incorporating local pairwise alignment information. Highly accurate, but slower.
- [MAFFT documentation](https://mafft.cbrc.jp/alignment/software/)

##### [`WITCH`](../modules/nf-core-modified/witch.nf):

- Parameters specified in [`conf/modules.config`](../conf/modules.config).
- See [WITCH documentation](https://github.com/c5shen/WITCH) for detailed description of options.
- **Known limitation (non-determinism):** WITCH selects its alignment backbone with an unseeded random sample, and the upstream `witch-msa` package exposes no option to set a random seed. As a result, WITCH alignments (and anything downstream of them) may vary slightly between runs on the same input. The proper long-term fix is an upstream `--seed` feature ([c5shen/WITCH](https://github.com/c5shen/WITCH)); we deliberately do not patch the third-party package internals.

##### [`FAMSA`](../modules/local/famsa.nf):

- Fast and accurate multiple sequence alignment algorithm optimized for large protein families
- Parameters specified in [`conf/modules.config`](../conf/modules.config)
- Default: Uses standard FAMSA parameters with automatic thread detection
- Common options (can be added to modules.config):
  - `-medoidtree`: Use medoid tree heuristic for faster alignment of very large families
  - `-gt upgma/sl/parttree`: Guide tree method selection (default: upgma)
  - `-t <n>`: Number of threads (automatically set from task.cpus)
- [FAMSA documentation](https://github.com/refresh-bio/FAMSA)

#### 6. `TRIM_SEQS`

##### [`CLIPKIT`](../modules/local/clipkit.nf):

- Defaults used. Custom parameters should be specified in [`conf/modules.config`](../conf/modules.config).
- [ClipKIT documentation](https://jlsteenwyk.com/ClipKIT/)

##### [`CIALIGN`](../modules/local/cialign.nf):

- Custom parameters specified in [`conf/modules.config`](../conf/modules.config).
- See [CIALIGN documentation](https://github.com/KatyBrown/CIAlign) for detailed description of options.
- `--crop_divergent_min_prop_ident=0.25 --crop_divergent_min_prop_nongap=0.25 --crop_ends --remove_insertions --insertion_min_size=5 --insertion_max_size=200 --remove_divergent --remove_divergent_minperc=0.15`

#### 7. `INFER_TREES`

##### [`FASTTREE`](../modules/nf-core-modified/fasttree.nf):

- Custom parameters specified in [`conf/modules.config`](../conf/modules.config).
- See [FastTree2 documentation](http://www.microbesonline.org/fasttree/) for detailed description of options.
- `-lg -gamma -bionj -pseudo -spr 4 -mlacc 2 -slownni`

##### [`IQTREE`](../modules/nf-core-modified/iqtree.nf):

- `tree_model`: Model of amino acid substitution to be used for phylogenetic inference. Specified in parameter-file.
- All other custom parameters should be specified in [`conf/modules.config`](../conf/modules.config).
- [IQ-TREE documentation](http://www.iqtree.org/)

#### 8. MiniNJ initialization

MiniNJ is run inside [`SPECIESRAX`](../modules/local/speciesrax.nf) as an
initialization-only pass. It uses all selected core gene trees under MPI, with
species-tree search and reconciliation disabled. NovelTree then replaces its
arbitrary serialized root using the trusted chronogram or explicit outgroups.

#### 9. [`SPECIESRAX`](../modules/local/speciesrax.nf):

##### **PLEASE** read the [SpeciesRax documentation](https://github.com/BenoitMorel/GeneRax/wiki/GeneRax) to GeneRax and SpeciesRax for a more detailed explanation, both of these options as well as other possible parameter specifications.

- The following parameters are specified within the [SpeciesRax module file](../modules/local/speciesrax.nf)
- `--strategy SKIP --si-estimate-bl`

- The following parameters are specified in [`conf/modules.config`](../conf/modules.config).
- `--rec-model UndatedDL --si-strategy SKIP --si-quartet-support`

Before constructing the SpeciesRax family file, NovelTree currently applies a
late-stage scaling guard to the final validated gene trees. Families must retain
at least 50% species occupancy, have no more than four mean copies among
represented species, have no more than eight copies in any one species, and
contain no more than four times the dataset species count in total leaves. This
retains multicopy information for duplication/loss-aware inference while
preventing exceptionally expanded families from dominating SpeciesRax runtime.
`speciesrax_family_selection.tsv` records every decision and exclusion reason,
and `speciesrax_selected_species_coverage.tsv` reports retained coverage by
species.

This is intentionally implemented within `SPECIESRAX` for now so existing
alignments and gene trees remain reusable. After the cutoffs have been validated
across production datasets, this classification will be incorporated into the
upstream orthogroup filtering and routing logic.

#### 10. [`GENERAX_PER_FAMILY`](../modules/local/generax_per_family.nf):

- The following parameters are specified within the [GeneRax per-family module file](../modules/local/generax_per_family.nf)
- `--prune-species-tree --reconciliation-samples 100`

- The following parameters are specified in [`conf/modules.config`](../conf/modules.config).
- `--rec-model UndatedDL --strategy SPR`

- [GeneRax documentation](https://github.com/BenoitMorel/GeneRax/wiki/GeneRax)

#### 11. [`GENERAX_PER_SPECIES`](../modules/local/generax_per_species.nf):

- The following parameters are specified within the [GeneRax per-species module file](../modules/local/generax_per_species.nf)
- `--prune-species-tree --reconciliation-samples 100 --per-species-rates`

- The following parameters are specified in [`conf/modules.config`](../conf/modules.config).
- `--rec-model UndatedDL --strategy SPR`

- [GeneRax documentation](https://github.com/BenoitMorel/GeneRax/wiki/GeneRax)

#### 12. [`PARSE_PHYLOHOGS`](../modules/local/parse_phylohogs.nf):

- Extracts ortholog/paralog pairs and hierarchical orthogroup (HOG) membership from GeneRax reconciliation output
- Inputs: GeneRax `_reconciliated.nhx` file + GeneRax-labeled species tree
- Outputs: `{OG}_orthologs.tsv`, `{OG}_paralogs.tsv`, `{OG}_hog_membership.tsv`, `spp_tree_node_lookup.tsv`

#### 13. [`PHYLO_PROFILES`](../modules/local/phylo_profiles.nf):

- Generates phylogenetic profiles from GeneRax per-species reconciliation outputs
- Summarizes gene duplication, loss, and speciation events across species
- No parameters required

#### 14a. [`BUILD_REFERENCE_CHRONOGRAM`](../modules/local/build_reference_chronogram.nf) _(zoogle mode only)_:

- Automatically builds a reference chronogram by querying TimeTree.org for pairwise divergence times
- Queries species-level taxids first, with genus-level fallback for missing pairs
- Uses greedy-drop to remove species with the most missing pairwise times until a complete distance matrix is obtained (requires ≥4 species)
- Constructs UPGMA tree via scipy, producing an ultrametric Newick tree
- **Requires**: `ncbi_email` parameter for NCBI Entrez queries
- Skipped when `reference_time_tree` is provided

#### 14b. [`TIME_CALIBRATE_SPECIES_TREE`](../modules/local/time_calibrate_species_tree.nf) _(zoogle mode only)_:

- Time-calibrates the inferred species tree against the reference chronogram (auto-built or user-provided)
- Matches shared taxa between the inferred and reference trees, then uses treePL penalized likelihood to date the species tree
- Applies `age_bracket` (default ±20%) around reference ages as min/max calibration bounds
- **Optional**: `reference_time_tree` parameter with path to reference timetree (Newick format). If not provided, a chronogram is auto-built from TimeTree.org.
- `time_calibration_method`: Method for calibration - "treePL" (default)

#### 15. [`DATE_GENE_FAMILY_TREES`](../modules/local/date_gene_family_trees.nf) _(zoogle mode only)_:

- Time-calibrates gene family trees using speciation node ages from the dated species tree
- Uses GeneRax reconciliation output (`_events.newick`) to identify speciation nodes — only speciation events are used as calibration points (duplications are excluded)
- Internal node labels (S/D) are stripped before dating to avoid confusing treePL
- Uses treePL with fixed `smooth=10` (no cross-validation pass) to avoid CV instability on large gene trees
- Calibration ages use `age_bracket` (default ±20%) around species tree node ages

#### 16. [`PROTEIN_PROPERTIES`](../modules/local/protein_properties.nf) _(zoogle mode only)_:

- Calculates amino acid composition and physicochemical properties for all gene families
- Computes 20 amino acid frequencies and properties (molecular weight, aromaticity, GRAVY, isoelectric point, etc.)
- No parameters required

#### 17. [`ZOOGLE`](../modules/local/zoogle.nf) _(zoogle mode only)_:

- Calculates phylogenetically-corrected protein distances using Mahalanobis distances
- **Universal centroid analysis**: For all gene families, computes Mahalanobis distance from each protein to the family centroid, with protein-level rank p-values and species-level permutation p-values
- **Reference analysis** _(optional)_: When `ref_species` is set (not `'none'`), performs permutation tests to identify proteins exceptionally (dis)similar to the reference species
- Processes all gene families with ≥4 proteins from ≥2 species (each with ≥2 proteins)
- Set `--ref_species none` for centroid-only analysis without a reference species
