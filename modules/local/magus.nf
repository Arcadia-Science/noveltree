process MAGUS {
    tag "$fasta"
    label 'process_mafft'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/vft-magus_0.1.0:0.0.1' :
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/vft-magus_0.1.0:0.0.1' :
        '' }"
    // TODO: address this issue (permission related errors) in future release
    containerOptions = "--user root"

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
    
    # Hacky fix to prevent segfault of VeryFastTree when running on small datasets
    nseqs=\$(grep ">" \$og | wc -l)
    if [[ \${nseqs} -gt 14 ]]
    then
        nthreads=${task.cpus}
    else
        nthreads=1
    fi
    
    magus \\
        -i $fasta \\
        -o \${prefix}_magus.fa \\
        --numprocs ${task.cpus} \\
        --numprocs_vft \${nthreads}
        -t veryfasttree
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magus-vft: v\$(grep "version=" MAGUS-VFT/setup.py | cut -f2 -d'"')
    END_VERSIONS
    """
}
