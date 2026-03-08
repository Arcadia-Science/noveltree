process PREQUAL {
    tag "${meta.og}"
    label 'process_single_nomem'

    // PREQUAL pairHMM is O(n_seq × max_len²). Estimate memory from family metadata.
    memory {
        def n = (meta?.n_seq ?: 10) as long
        def L = (meta?.max_len ?: 500) as long
        // Rough estimate: 2 × L² × 8 bytes per sequence pair, plus 1 GB overhead
        def estimated_gb = Math.max(2L, (long)(2L * n * L * L * 8 / (1024 * 1024 * 1024)) + 1L)
        def capped_gb = Math.min(estimated_gb, 64L)
        check_max((capped_gb as int).GB * task.attempt, 'memory')
    }

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
