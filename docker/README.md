# Arcadia-Science/noveltree Docker images

This folder includes the Dockerfiles used for the custom modules of the pipeline. The final versions of all containers are hosted on [Arcadia Science's DockerHub account](https://hub.docker.com/u/arcadiascience).

## Base images

The Dockerfiles we include in this folder use three separate base images as starting points: [ubuntu:20.04](https://hub.docker.com/layers/library/ubuntu/20.04/images/sha256-3246518d9735254519e1b2ff35f95686e4a5011c90c85344c1f38df7bae9dd37?context=explore), [pyton:3.9-slim](https://hub.docker.com/layers/library/python/3.9-slim/images/sha256-b370e60efdfcc5fcb0a080c0905bbcbeb1060db3ce07c3ea0e830b0d4a17f758) and [arcadiascience/rbase_4.2.2](https://hub.docker.com/r/arcadiascience/rbase_4.2.2). The final image is a custom built image by us that comes with R v4.2.2.

## Naming conventions

The Docker image names follow this convention: `<BASE_NAME>:<SEMANTIC_VERSIONING>`, all in lower-case.

`<BASE_NAME>` part follows these guides:
* Describe what the image essentially does or represents. If it's a web service with Node.js and Nginx, it might be something like `node-nginx`. If it only includes `Asteroid`, it might be something like `asteroid`. If an image uses way too many libraries, a descriptive name is preferred (i.e. `select_mcl_inflation_params`).
  * If multiple packages need to be represented, they should be separated by `-` (i.e. `node-nginx`)
* The main package's version should be represented in the image name. For example, if you're shipping Python 3.8, you could use python3.8 as part of the name or tag. If a version number is available, the first 8 characters of the commit SHA can be used. If none of that is available, you can use the date stamp. For instance:
  * If an image is using ClipKIT v2.1.1, the base name would be `clipkit_2.1.1`
  * If an image is using Asteroid (only available through GitHub) at commit SHA `3aae117df3353c28c6a07d58a4c8b0ab290f974f`, the base name would be `asteroid_3aae117d`

`<SEMANTIC_VERSIONING>` is a widely accepted standard. A version like `1.2.3` stands for `MAJOR.MINOR.PATCH``. This provides clear semantics on the kind of changes that happened. This version is different from the underlying package version and it reflects the level of stability of the image from our perspective.

## Module images

| Process name | Docker image | Container |
| --- | --- | --- |
| `SAMPLESHEET_CHECK` | [python 3.9](./python/) | `arcadiascience/python_3.9` |
| `DOWNLOAD_INPUT`, `TRANSDECODER`, `FILTER_ISOFORMS`, `PREPROCESS_PROTEOME` | [preprocess_proteomes v1.1.0](./preprocess_proteomes/) | `arcadiascience/preprocess_proteomes:1.1.0` |
| `RENAME_FASTAS` | [R v4.2.2](./rbase/) | `arcadiascience/rbase_4.2.2:1.0.0` |
| `ANNOTATE_UNIPROT` | [bioservices v1.10.0](./bioservices/) | `arcadiascience/bioservices_1.10.0:1.0.0` |
| `ASTEROID` | [Asteroid @ commit 3aae117d](./asteroid/) | `arcadiascience/asteroid_3aae117d-disco_20e10c33:1.0.0` |
| `CLIPKIT` | [ClipKIT v2.1.1](./clipkit/) | `arcadiascience/clipkit_2.1.1-seqmagick_0.8.4:1.0.0` |
| `CIALIGN` | [CIAlign v1.1.0](./cialign/) | `arcadiascience/cialign_1.1.0:1.0.0` |
| `COGEQC` | [cogeqc v1.2.1](./cogeqc/) | `arcadiascience/cogeqc_1.2.1:1.0.0` |
| `FAMSA` | [FAMSA 2.5.0 @ commit 5a326d5 (container name `famsa_2.0.0` is historical)](./famsa/) | `arcadiascience/famsa_2.0.0:1.0.0` |
| `FASTTREE` | [FastTree v2.1.11](./fasttree/) | `arcadiascience/fasttree_2.1.11:1.0.0` |
| `GENERAX_PER_FAMILY`, `GENERAX_PER_SPECIES` | [GeneRax @ commit 56f3ed0](./generax/) | `arcadiascience/generax_56f3ed0:1.1.3` |
| `SPECIESRAX` | [GeneRax @ commit 56f3ed0 + Biopython 1.83](./generax/) | `arcadiascience/generax_56f3ed0:1.1.4` |
| `IQTREE` | [IQ-TREE v2.2.0.5](./iqtree/) | `arcadiascience/iqtree_2.2.0.5:1.0.0` |
| `ORTHOFINDER_PREP`, `ORTHOFINDER_MCL` | [OrthoFinder v2.5.4](./orthofinder/) | `arcadiascience/orthofinder_2.5.4:1.0.0` |
| `PARSE_PHYLOHOGS`, `PHYLO_PROFILES` | [phylo_profiles v1.0.0](./phylo_profiles/) | `arcadiascience/phylo_profiles:1.0.0` |
| `PROTEIN_PROPERTIES` | [protein_properties v1.0.0](./physicochemical_props/) | `arcadiascience/protein_properties:1.0.0` |
| `SELECT_INFLATION` | [select_mcl_inflation_params](./select_mcl_inflation_params/) | `arcadiascience/select_mcl_inflation_params_08302023:1.0.0` |
| `WITCH` | [WITCH v1.0.10](./witch/) | `arcadiascience/witch_1.0.10:1.0.0` |
| `BUILD_REFERENCE_CHRONOGRAM` | [build_reference_chronogram v1.0.0](./build_reference_chronogram/) | `arcadiascience/build_reference_chronogram:1.0.0` |
| `ZOOGLE_ANALYSIS`, `TIME_CALIBRATE_SPECIES_TREE`, `DATE_GENE_FAMILY_TREES` | [zoogle v1.2.0](./zoogle/) | `arcadiascience/zoogle:1.2.0` |

## Reproducibility and upstream availability

Several Dockerfiles build tools directly from third-party GitHub repositories (e.g. Asteroid, DISCO, GeneRax, CIAlign, IQ-TREE 2, FAMSA). To keep builds reproducible, each of these is pinned to a specific upstream commit SHA or release tag rather than a moving branch. There are two independent layers of protection:

1. **Pinned build inputs** — the exact commit/tag in each Dockerfile means a rebuild produces the same tool version as the published image.
2. **Pre-built images on DockerHub** — the authoritative, ready-to-run artifacts already live on [Arcadia Science's DockerHub account](https://hub.docker.com/u/arcadiascience) and are pulled automatically by the pipeline. These do not depend on the upstream repositories remaining online.

If an upstream repository is ever renamed, deleted, or made unavailable, **the published DockerHub image remains the fallback** — the pipeline continues to run without any change. Rebuilding that image from scratch would then require vendoring the tool's source (e.g. from a fork or an archived copy) and repointing the Dockerfile at it. We do not currently maintain Arcadia forks/mirrors of these upstream repositories.
