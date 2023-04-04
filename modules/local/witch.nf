process WITCH {
    tag "$fasta"
    label 'process_magus'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/witch_0.3.0:0.0.1' :
        '' }"
    // TODO: address this issue (permission related errors) in future release
    containerOptions = "--user root"
    
    stageInMode = 'copy'
    publishDir(
        path: "${params.outdir}/witch_alignments",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file(fasta)

    output:
    path("**_witch.fa")          , emit: msas
    path("**_witch_cleaned.fa")  , emit: cleaned_msas
    path("*")                    , emit: results
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def min_len = params.min_ungapped_length ?: '20'
    """
    prefix=\$(basename "${fasta}" .fa)

    # WITCH does unfortunately need a little assistance for smaller alignments,
    # where including too many sequences within the backbone alignment will 
    # lead to situations where no sequences will be within 25% of the median
    # length
    
    # So, to do this, we'll use a slightly modified version of the backbone.py script
    # that will automatically rescale the minimum and maximum length of sequences 
    # that are considered "full length"
    #ntax=\$(grep ">" ${fasta} | wc -l)

    # Be sure to remove any non-standard amino acid codes in the input sequences, as this 
    # can cause errors downstream and in parsing. 
    sed -E -i '/>/!s/U//g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O//g' ${fasta} # pyrrolysine
    
    #set +e # Turn off error recognition
    python3 /WITCH/witch.py \\
        -i ${fasta} \\
        -d alignments \\
        -t ${task.cpus} \\
        --graphtraceoptimize true \\
        --molecule amino \\
        $args
    
    # By default, witch only masks (and removes) singleton columns - we still would 
    # like to make sure that we are removing sequences with too many positions removed. 
    # So, here we're using awk to remove sequences with fewer than params.min_ungapped_length
    # AA remaining once masked. 
    awk -v N=${min_len} -F "" \
        '{ s=0; for (i=1; i<=NF; i++) if (\$i != "-") s++ } /^>/ { if (s >= N || NR == 1) print; \
        if (s >= N) f=1; else f=0 } !/^>/ { if (f) print }' \
        alignments/merged.fasta.masked > final_masked.fasta
        
    # Reorganize results for publishing
    mkdir original_alignments
    mkdir cleaned_alignments
    mv alignments/merged.fasta original_alignments/\${prefix}_witch.fa
    mv final_masked.fasta cleaned_alignments/\${prefix}_witch_cleaned.fa
    rm -r alignments/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        witch: v\$(python3 /WITCH/witch.py -v | cut -f2 -d " ")
    END_VERSIONS
    """
}
