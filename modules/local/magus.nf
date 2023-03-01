process MAGUS {
    tag "$fasta"
    label 'process_mafft'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/vft-magus_0.1.0:0.0.1' :
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/vft-magus_0.1.0:0.0.1' :
        '' }"

    input:
    path(fasta)
    val num_spp_tree_msas // Purely utilitarian: used to prioritize gene families used in species tree inference.

    output:
    path("*_magus.fa")  , emit: msas
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prefix=\$(basename "${fasta}" .fa)
    magus \\
        -i $fasta \\
        -o \${prefix}_magus.fa \\
        --numprocs ${task.cpus} \\
        -t veryfasttree
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magus-vft: 0.1.0
    END_VERSIONS
    """
}
