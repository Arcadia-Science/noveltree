process PREQUAL {
    tag "${meta.og}"
    label 'process_single'

    container 'arcadiascience/prequal:1.0.0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.og}_prequal.fa"), emit: masked
    path "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prequal "${fasta}" ${args}

    # PREQUAL produces input.fasta.filtered; rename to expected output
    mv "${fasta}.filtered" "${meta.og}_prequal.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prequal: \$(prequal 2>&1 | head -1 | grep -oP '\\d+\\.\\d+' || echo "1.02")
    END_VERSIONS
    """
}
