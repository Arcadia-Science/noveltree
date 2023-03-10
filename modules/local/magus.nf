process MAGUS {
    tag "$fasta"
    label 'process_magus'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/magus_0.1.0:0.0.1' :
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/magus_0.1.0:0.0.1' :
        '' }"
    // TODO: address this issue (permission related errors) in future release
    containerOptions = "--user root"

    publishDir(
        path: "${params.outdir}/magus_alignments",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
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
        -i ${fasta} \\
        -o \${prefix}_magus.fa \\
        --numprocs ${task.cpus} \\
        $args
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magus: v\$(grep "version=" /MAGUS/setup.py | cut -f2 -d'"')
    END_VERSIONS
    """
}