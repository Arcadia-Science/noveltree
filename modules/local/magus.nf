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
    
    if [[ \$ntax -ge 200 ]]; then
        # Skeleton size = 20% of taxa
        skelsize=\$(printf "%.0f" \$(echo "scale=2;\$ntax*0.20" | bc))"
        decompskeletonsize="--decompskeletonsize \$skelsize"
        
        # Max subset size = 25% of skeleton size
        subsetsize=\$(printf "%.0f" \$(echo "scale=2;\$skelsize*0.25" | bc))
        maxsubsetsize="--maxsubsetsize \$subsetsize"
        
        # Max number of subsets = 20 - max subset size will typically be 'favored'
        maxnumsubsets="--maxnumsubsets 20"
        
        # Mafft backbone alignment size = 20% of taxa
        mafftsize="-m \$(printf "%.0f" \$(echo "scale=2;\$subsetsize*0.20" | bc))"
        
        # Number of Mafft runs = 5
        mafftruns="-r 5"
        
        graphbuildhmmextend="--graphbuildhmmextend true"
    elif [[ \$ntax -le 199 && \$ntax -ge 100 ]]; then
        skelsize=\$(printf "%.0f" \$(echo "scale=2;\$ntax*0.30" | bc))
        decompskeletonsize="--decompskeletonsize \$skelsize"
        subsetsize=\$(printf "%.0f" \$(echo "scale=2;\$skelsize*0.25" | bc))
        maxsubsetsize="--maxsubsetsize \$subsetsize"
        maxnumsubsets="--maxnumsubsets 20"
        mafftsize="-m \$(printf "%.0f" \$(echo "scale=2;\$subsetsize*0.20" | bc))"
        mafftruns="-r 5"
        graphbuildhmmextend="--graphbuildhmmextend true"
    elif [[ \$ntax -le 99 && \$ntax -ge 20 ]]; then
        skelsize=\$(printf "%.0f" $(echo "scale=2;\$ntax*0.50" | bc))
        decompskeletonsize="--decompskeletonsize \$skelsize"
        subsetsize=\$(printf "%.0f" \$(echo "scale=2;\$skelsize*0.40" | bc))
        maxsubsetsize="--maxsubsetsize \$subsetsize"
        maxnumsubsets="--maxnumsubsets 10"
        mafftsize="-m \$(printf "%.0f" \$(echo "scale=2;\$subsetsize*0.70" | bc))"
        mafftruns="-r 3"
        graphbuildhmmextend="--graphbuildhmmextend true"
    elif [[ \$ntax -le 19 && \$ntax -ge 10 ]]; then
        skelsize=$(printf "%.0f" $(echo "scale=2;\$ntax*0.65" | bc))
        decompskeletonsize="--decompskeletonsize \$skelsize"
        subsetsize=\$(printf "%.0f" \$(echo "scale=2;\$skelsize*0.50" | bc))
        maxsubsetsize="--maxsubsetsize \$subsetsize"
        maxnumsubsets="--maxnumsubsets 5"
        mafftsize="-m \$(printf "%.0f" \$(echo "scale=2;\$subsetsize*0.70" | bc))"
        mafftruns="-r 3"
        graphbuildhmmextend="--graphbuildhmmextend true"
    else
        skelsize=\$ntax
        decompskeletonsize="--decompskeletonsize \$skelsize"
        maxsubsetsize="--maxsubsetsize 3"
        maxnumsubsets="--maxnumsubsets 3"
        mafftsize="-m 3"
        mafftruns="-r 3"
        graphbuildhmmextend="--graphbuildhmmextend false"
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
