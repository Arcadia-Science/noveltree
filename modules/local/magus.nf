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
    
    # Prevent needless excess subsetting by dynamically specifying here
    ntax=\$(grep ">" ${fasta} | wc -l)
    if [ \$ntax -ge 100 ]; then
        decompskeletonsize="--decompskeletonsize 101"
        maxsubsetsize="--maxsubsetsize 50"
        maxnumsubsets="--maxnumsubsets 20"
        graphbuildhmmextend="--graphbuildhmmextend true"
        mafftsize="-m 25"
    elif [ \$ntax -le 100 && \$ntax -ge 25 ]; then 
        decompskeletonsize="--decompskeletonsize 25"
        maxsubsetsize="--maxsubsetsize 10"
        maxnumsubsets="--maxnumsubsets 10"
        graphbuildhmmextend="--graphbuildhmmextend true"
        mafftsize="-m 10"
    else [ \$ntax -le 24 && \$ntax -ge 10 ]
        decompskeletonsize="--decompskeletonsize 24"
        maxsubsetsize="--maxsubsetsize 25"
        maxnumsubsets="--maxnumsubsets 1"
        graphbuildhmmextend="--graphbuildhmmextend false"
        mafftsize="-m \${ntax}"
    fi
        
    magus \\
        -i ${fasta} \\
        -o \${prefix}_magus.fa \\
        --numprocs ${task.cpus} \\
        \${decompskeletonsize} \\
        \${maxsubsetsize} \\
        \${maxnumsubsets} \\
        \${graphbuildhmmextend} \\
        \${mafftsize} 
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magus: v\$(grep "version=" /MAGUS/setup.py | cut -f2 -d'"')
    END_VERSIONS
    """
}
