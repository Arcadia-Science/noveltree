process BUNDLE_BLAST_RESULTS {
    tag "$bundle_id"
    label 'process_single'

    storeDir "${params.outdir}/blast_bundles"

    container 'arcadiascience/orthofinder_2.5.4:1.0.0'

    input:
    tuple val(bundle_id), path(blast_files)

    output:
    path "BlastBundle_${bundle_id}.tar", emit: archive

    script:
    """
    # DIAMOND output is already gzip-compressed. Store it in an uncompressed
    # tar archive so that large datasets are staged into OrthoFinder as one
    # archive per query species instead of one input per species pair.
    tar -cf BlastBundle_${bundle_id}.tar *.txt.gz
    """
}
