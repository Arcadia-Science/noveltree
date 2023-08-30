# Arcadia-Science/noveltree Docker images

This folder includes the Dockerfiles used for the custome modules of the pipeline. The final versions of all containers are hosted on [Arcadia Science's DockerHub account](https://hub.docker.com/u/arcadiascience).

## Base images

The Dockerfiles we include in this folder use 2 separate base images as starting points: [ubuntu:20.04](https://hub.docker.com/layers/library/ubuntu/20.04/images/sha256-3246518d9735254519e1b2ff35f95686e4a5011c90c85344c1f38df7bae9dd37?context=explore) and [pyton:3.9-slim](https://hub.docker.com/layers/library/python/3.9-slim/images/sha256-b370e60efdfcc5fcb0a080c0905bbcbeb1060db3ce07c3ea0e830b0d4a17f758).

## Naming conventions

The Docker image names follow this convention: `<BASE_NAME>:<SEMANTIC_VERSIONING>`, all in lower-case.

`<BASE_NAME>` part follows these guides:
* Describe what the image essentially does or represents. If it's a web service with Node.js and Nginx, it might be something like `node-nginx`. If it only includes `Asteroid`, it might be something like `asteroid`. If an image uses way too many libraries, a descriptive name is preferred (i.e. `slect_mcl_inflation_params`).
  * If multiple packages need to be represented, they should be separated by `-` (i.e. `node-nginx`)
* The main package's version should be represented in the image name. For example, if you're shipping Python 3.8, you could use python3.8 as part of the name or tag. If a version number is available, the first 8 characters of the commit SHA can be used. If none of that is available, you can use the date stamp. For instance:
  * If an image is using Clipkit v1.3.0, the base name would be `clipkit_1.3.0`
  * If an image is using Asteroid (only available through GitHub) at commit SHA `3aae117df3353c28c6a07d58a4c8b0ab290f974f`, the base name would be `asteroid_3aae117d`

`<SEMANTIC_VERSIONING>` is a widely accepted standard. A version like `1.2.3` stands for `MAJOR.MINOR.PATCH``. This provides clear semantics on the kind of changes that happened. This version is different from the underlying package version and it reflects the level of stability of the image from our perspective.

## Module images

| Process name                         | Docker image                                                         | Hosted image location                                                                                                    |
| ------------------------------------ | -------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `ASTEROID`                           | [Asteroid @ commit 3aae117d](./asteroid/)                             | [arcadiascience/asteroid_3aae117d-disco_20e10c33](https://hub.docker.com/r/arcadiascience/asteroid_3aae117d-disco_20e10c33)
| `ANNOTATE_UNIPROT`                   | [bioservices v1.10.0](./bioservices/)                             | [arcadiascience/bioservices_1.10.0](https://hub.docker.com/r/arcadiascience/bioservices_1.10.0)                       |
| `CIALIGN`                            | [CIAlign v1.1.0](./cialign/)                                         | [arcadiascience/cialign_1.1.0](https://hub.docker.com/r/arcadiascience/cialign_1.1.0)                             |
| `CLIPKIT`                            | [clipkit v1.3.0](./clipkit/)                                         | [arcadiascience/clipkit_1.3.0-seqmagick_0.8.4](https://hub.docker.com/r/arcadiascience/clipkit_1.3.0-seqmagick_0.8.4)                             |
| `COGEQC`                             | [cogeqc v1.2.1](./cogeqc/)                                           | [arcadiascience/cogeqc_1.2.1](https://hub.docker.com/r/arcadiascience/cogeqc_1.2.1)                               |
| `FASTTREE`                 | [FastTree v2.1.11](./fasttree/)                               | [arcadiascience/fasttree_2.1.11](https://hub.docker.com/r/arcadiascience/fasttree_2.1.11)                         |
| `GENERAX/SPECIESRAX`                 | [GeneRax @ commit 19604b71](./generax/)                               | [arcadiascience/generax_19604b71](https://hub.docker.com/r/arcadiascience/generax_19604b71)                         |
| `IQTREE`                             | [IQ-TREE v2.2.0.5](./iqtree/)                                        | [arcadiascience/iqtree_2.2.0.5](https://hub.docker.com/r/arcadiascience/iqtree_2.2.0.5)                        |
| `ORTHOFINDER` (prep, mcl, phylohogs) | [OrthoFinder v2.5.4](./orthofinder)                                  | [arcadiascience/orthofinder_2.5.4](https://hub.docker.com/r/arcadiascience/orthofinder_2.5.4)                     |
| `FILTER_ORTHOGROUPS`                 | [R v4.2.2](./rbase/)                                                 | [arcadiascience/rbase_4.2.2](https://hub.docker.com/r/arcadiascience/rbase_4.2.2)                                 |
| `SELECT_INFLATION`                   | [select_mcl_inflation_params created on 08/30/2023](./select_mcl_inflation_params/) | [arcadiascience/select_mcl_inflation_params_08302023](https://hub.docker.com/r/arcadiascience/select_mcl_inflation_params_08302023) |
| `WITCH`                   | [witch v0.3.0](./witch/) | [arcadiascience/witch_0.3.0](https://hub.docker.com/r/arcadiascience/witch_0.3.0) |
