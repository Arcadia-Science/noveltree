process TRANSDECODER {
    tag "${meta.id}"
    label 'process_medium'

    container 'arcadiascience/preprocess_proteomes:1.0.0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_translated.fasta"), emit: translated
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Decompress if gzipped
    input_file="${fasta}"
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > input.fasta
        input_file="input.fasta"
    fi

    TransDecoder.LongOrfs -t "\${input_file}"
    TransDecoder.Predict -t "\${input_file}" --single_best_only

    cp "\${input_file}.transdecoder.pep" "${meta.id}_translated.fasta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(TransDecoder.LongOrfs --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo "5.7.1")
    END_VERSIONS
    """
}
