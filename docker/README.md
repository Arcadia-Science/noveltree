# Arcadia-Science/phylorthology Docker images

This folder includes the Dockerfile used for the custome modules of the pipeline. The final versions of all containers are hosted on [Arcadia Science's DockerHub account](https://hub.docker.com/u/arcadiascience).

## Images

| Process name                         | Docker image                                                         | Image location                                                                                                    |
| ------------------------------------ | -------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `ASTEROID`                           | [Asteroid @ commit 3aae117](./asteroid/)                             | [arcadiascience/asteroid_3aae117](https://hub.docker.com/r/arcadiascience/asteroid_3aae117)                       |
| `CIALIGN`                            | [CIAlign v1.1.0](./cialign/)                                         | [arcadiascience/clipkit_1.3.0](https://hub.docker.com/r/arcadiascience/clipkit_1.3.0)                             |
| `CLIPKIT`                            | [clipkit v1.3.0](./clipkit/)                                         | [arcadiascience/clipkit_1.3.0](https://hub.docker.com/r/arcadiascience/clipkit_1.3.0)                             |
| `COGEQC`                             | [cogeqc v1.2.1](./cogeqc/)                                           | [arcadiascience/cogeqc_1.2.1](https://hub.docker.com/r/arcadiascience/cogeqc_1.2.1)                               |
| `GENERAX/SPECIESRAX`                 | [GeneRax @ commit 19604b7](./generax/)                               | [arcadiascience/generax_19604b7](https://hub.docker.com/r/arcadiascience/generax_19604b7)                         |
| `IQTREE`                             | [IQ-TREE v2.2.0.5](./iqtree/)                                        | [arcadiascience/iqtree\_:_2.2.0.5](https://hub.docker.com/r/arcadiascience/iqtree_2.2.0.5)                        |
| `ORTHOFINDER` (prep, mcl, phylohogs) | [OrthoFinder v2.5.4](./orthofinder)                                  | [arcadiascience/orthofinder_2.5.4](https://hub.docker.com/r/arcadiascience/orthofinder_2.5.4)                     |
| `FILTER_ORTHOGROUPS`                 | [R v4.2.2](./rbase/)                                                 | [arcadiascience/rbase_4.2.2](https://hub.docker.com/r/arcadiascience/rbase_4.2.2)                                 |
| `SELECT_INFLATION`                   | [select_mcl_inflation_params v0.0.1](./select_mcl_inflation_params/) | [arcadiascience/select_mcl_inflation_params](https://hub.docker.com/r/arcadiascience/select_mcl_inflation_params) |
| `ANNOTATE_UNIPROT`                   | [UniProt.ws v2.38.1](./uniprotws/)                                   | [arcadiascience/uniprotws_2.38.1](https://hub.docker.com/r/arcadiascience/uniprotws_2.38.1)                       |
