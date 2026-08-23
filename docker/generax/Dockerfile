# Preserve the exact GeneRax binary already used by the production workflow.
# This SpeciesRax-only derivative adds the Newick tooling needed to transfer a
# trusted reference root without rebuilding or changing GeneRax itself.
FROM arcadiascience/generax_56f3ed0:1.1.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends python3-pip && \
    rm -rf /var/lib/apt/lists/*

# SpeciesRax computes MiniNJ and transfers a trusted reference root within one
# task; Biopython provides the validated Newick rooting implementation.
RUN python3 -m pip install --no-cache-dir biopython==1.83
