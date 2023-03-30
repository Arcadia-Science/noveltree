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
    path(fasta)

    output:
    path("**_witch.fa")          , emit: msas
    path("**_witch_cleaned.fa")  , emit: cleaned_msas
    path("*")                    , emit: results
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prefix=\$(basename "${fasta}" .fa)

    # WITCH does unfortunately need a little assistance for smaller alignments,
    # where including too many sequences within the backbone alignment will 
    # lead to situations where no sequences will be within 25% of the median
    # length
    
    # So, we'll try to nudge this along with as little messing about as possible
    # We need to change the config file to specify the backbone size and the 
    # backbone size (at the time there is not a way to do so via a commandline argument)
    ntax=\$(grep ">" ${fasta} | wc -l)

    # Be sure to remove any non-standard amino acid codes in the input sequences, as this 
    # can cause errors downstream and in parsing. 
    sed -E -i '/>/!s/U/-/g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O/-/g' ${fasta} # pyrrolysine
    
    set +e # Turn off error recognition
    python3 /WITCH/witch.py \\
        -i ${fasta} \\
        -d alignments \\
        -t ${task.cpus} \\
        --graphtraceoptimize true \\
        --molecule amino \\
        $args
        
    set -e # Turn back on - if this dies due to the length threshold being too stringent then retry
    if [[ ! -s alignments/merged.fasta.masked ]]; then
        sed -i "s/backbone_threshold = 0.25/backbone_threshold = 0.75/g" /WITCH/gcmm/backbone.py
        size=\$(echo "\$ntax" | awk '{printf "%.0f", \$0*0.25}')
        skelsize=\$((size < 2 ? 2 : size))
        sed -i "s/backbone_size =/backbone_size = \${skelsize}/g" /WITCH/main.config
        python3 /WITCH/witch.py \\
            -i ${fasta} \\
            -d alignments \\
            -t ${task.cpus} \\
            --graphtraceoptimize true \\
            --molecule amino \\
            $args
    fi
    
    # Reorganize results for publishing
    mkdir original_alignments
    mkdir cleaned_alignments
    mv alignments/merged.fasta original_alignments/\${prefix}_witch.fa
    mv alignments/merged.fasta.masked cleaned_alignments/\${prefix}_witch_cleaned.fa
    rm -r alignments/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        witch: v\$(python3 /WITCH/witch.py -v | cut -f2 -d " ")
    END_VERSIONS
    """
}
