process MAFFT {
    tag "$meta.og"
    label 'process_highthread'

    conda (params.enable_conda ? "bioconda::mafft=7.490" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.490--h779adbc_0':
        'quay.io/biocontainers/mafft:7.490--h779adbc_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*-mafft.fa"), emit: msas
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.og}"
    """
    mafft \\
        --thread ${task.cpus} \\
        --localpair \\
        --maxiterate 1000 \\
        --anysymbol \\
        ${fasta} \\
        > ${prefix}-mafft.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
