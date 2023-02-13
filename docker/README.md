# Arcadia-Science/phylorthology Docker images

This folder includes the Dockerfile used for the custome modules of the pipeline. The final versions of all containers are hosted on [Arcadia Science's DockerHub account](https://hub.docker.com/u/arcadiascience).

## Images

| Process name  | Docker image  | Previously used |
| ------------- | ------------- | ------------- |
| `ASTEROID`  | [Asteroid @ commit 3aae117][]  | `austinhpatton123/asteroid:1.0.0` |
| `CLIPKIT`  | [clipkit v1.3.0][]  | `austinhpatton123/clipkit` |
| `COGEQC`  | [cogeqc v1.2.1][]  | `austinhpatton123/cogeqc-1.2.0_r-4.2.2` |
| `GENERAX/SPECIESRAX`  | [GeneRax @ commit 19604b7][]  | `quay.io/biocontainers/generax:2.0.4--h19e7193\_0` |
| `ORTHOFINDER` (prep, mcl, phylohogs)  | [OrthoFinder v2.5.4][]  | `quay.io/biocontainers/generax:2.0.4--h19e7193\_0` |
| `FILTER_ORTHOGROUPS`  | [R v4.2.2][]  | `austinhpatton123/cogeqc-1.2.0_r-4.2.2` |
| `SELECT_INFLATION`  | [select_mcl_inflation_params v0.0.1][]  | `austinhpatton123/select_mcl_inflation_r-4.2.2_elbow_tidy_reshape_cowplot` |
| `ANNOTATE_UNIPROT`  | [UniProt.ws v2.38.1][]  | `austinhpatton123/r-4.2.2_uniprot.ws:2.38.0` |
