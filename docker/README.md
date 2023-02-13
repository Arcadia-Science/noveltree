# Arcadia-Science/phylorthology Docker images

This folder includes the Dockerfile used for the custome modules of the pipeline. The final versions of all containers are hosted on [Arcadia Science's DockerHub account](https://hub.docker.com/u/arcadiascience).

## Images

| Process name  | Docker image  | Image location |
| ------------- | ------------- | ------------- |
| `ASTEROID`  | [Asteroid @ commit 3aae117](./asteroid/)  | `arcadiascience/asteroid:3aae117` |
| `CLIPKIT`  | [clipkit v1.3.0](./clipkit/)  | `arcadiascience/clipkit:1.3.0` |
| `COGEQC`  | [cogeqc v1.2.1](./cogeqc/)  | `arcadiascience/cogeqc:1.2.1` |
| `GENERAX/SPECIESRAX`  | [GeneRax @ commit 19604b7](./generax/)  | `arcadiascience/generax:19604b7` |
| `ORTHOFINDER` (prep, mcl, phylohogs)  | [OrthoFinder v2.5.4](./orthofinder) | `arcadiascience/orthofinder:2.5.4` |
| `FILTER_ORTHOGROUPS`  | [R v4.2.2](./rbase/)  | `arcadiascience/rbase:4.2.2` |
| `SELECT_INFLATION`  | [select_mcl_inflation_params v0.0.1](./select_mcl_inflation_params/)  | `arcadiascience/select_mcl_inflation_params:0.0.1` |
| `ANNOTATE_UNIPROT`  | [UniProt.ws v2.38.1](./uniprotws/)  | `arcadiascience/uniprotws:2.38.1` |
