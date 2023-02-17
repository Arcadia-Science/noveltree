process MAFFT {
    // Modified from nf-core to:
    // 1) remove unnecessary meta input
    // 2) use MAFFTs built in automatic thread scaling to improve memory efficiency
    // 3) wait for MSAs of species tree gene families to finish prior to running on remainder
    tag "$fasta"
    label 'process_mafft'

    conda (params.enable_conda ? "bioconda::mafft=7.490" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.490--h779adbc_0':
        'quay.io/biocontainers/mafft:7.490--h779adbc_0' }"

    input:
    path(fasta)
    val num_spp_tree_msas // Purely utilitarian: used to prioritize gene families used in species tree inference.

    output:
    path("*_mafft.fa")  , emit: msas
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prefix=\$(basename "${fasta}" .fa)
    mafft \\
        --thread ${task.cpus} \\
        ${args} ${fasta} > \${prefix}_mafft.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
