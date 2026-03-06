process FILTER_ISOFORMS {
    tag "${meta.id}"
    label 'process_single'

    container 'arcadiascience/preprocess_proteomes:1.0.0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_isofiltered.fasta"), emit: filtered
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_isoforms.py "${fasta}" "${meta.id}_isofiltered.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_isoforms: 1.0.0
    END_VERSIONS
    """
}
