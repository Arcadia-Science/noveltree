process BUNDLE_SPECIESRAX_INPUTS {
    tag "SpeciesRax input shard ${shard_id}"
    label 'process_bundle'

    container 'arcadiascience/orthofinder_2.5.4:1.0.0'

    input:
    tuple val(shard_id), path(map_links), path(gene_trees)

    output:
    path "speciesrax_inputs_${shard_id}.tar", emit: archive

    script:
    """
    bundle_speciesrax_inputs.py --input-dir . --output speciesrax_inputs_${shard_id}.tar --manifest-name speciesrax_inputs_${shard_id}.tsv
    """
}
