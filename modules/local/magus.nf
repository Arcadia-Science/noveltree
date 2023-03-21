process MAGUS {
    tag "$fasta"
    label 'process_magus'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/magus_0.1.0:0.0.1' :
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
    if [[ \$ntax -ge 101 ]]; then
        decompskeletonsize="--decompskeletonsize 101"
        maxsubsetsize="--maxsubsetsize 50"
        maxnumsubsets="--maxnumsubsets 20"
        graphbuildhmmextend="--graphbuildhmmextend true"
        mafftsize="-m 25"
        mafftruns="-r 10"
    elif [[ \$ntax -le 100 && \$ntax -ge 25 ]]; then
        decompskeletonsize="--decompskeletonsize 25"
        maxsubsetsize="--maxsubsetsize 10"
        maxnumsubsets="--maxnumsubsets 10"
        graphbuildhmmextend="--graphbuildhmmextend true"
        mafftsize="-m 10"
        mafftruns="-r 5"
    elif [[ \$ntax -le 24 && \$ntax -ge 10 ]]; then
        decompskeletonsize="--decompskeletonsize 15"
        maxsubsetsize="--maxsubsetsize 10"
        maxnumsubsets="--maxnumsubsets 2"
        graphbuildhmmextend="--graphbuildhmmextend false"
        mafftsize="-m 5"
        mafftruns="-r 2"
    else
        decompskeletonsize="--decompskeletonsize 9"
        maxsubsetsize="--maxsubsetsize 5"
        maxnumsubsets="--maxnumsubsets 2"
        graphbuildhmmextend="--graphbuildhmmextend false"
        mafftsize="-m 4"
        mafftruns="-r 1"
    fi

    magus \\
        -i ${fasta} \\
        -o \${prefix}_magus.fa \\
        --numprocs ${task.cpus} \\
        \${decompskeletonsize} \\
        \${maxsubsetsize} \\
        \${maxnumsubsets} \\
        \${graphbuildhmmextend} \\
        \${mafftsize} \\
        \${mafftruns}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magus: v\$(grep "version=" /MAGUS/setup.py | cut -f2 -d'"')
    END_VERSIONS
    """
}
