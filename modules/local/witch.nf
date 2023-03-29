process WITCH {
    tag "$fasta"
    label 'process_magus'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/witch_0.3.0:0.0.1' :
        '' }"
    // TODO: address this issue (permission related errors) in future release
    containerOptions = "--user root"

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

    python3 WITCH/witch.py \\
        -i \${prefix}.fa \\
        -d alignments \\
        --graphtraceoptimize true \\
        --molecule amino \\
        $args
    
    # Reorganize results for publishing
    mkdir original_alignments
    mkdir cleaned_alignments
    mv alignments/merged.fasta original_alignments/\${prefix}_witch.fa
    mv alignments/merged.fasta.masked cleaned_alignments/\${prefix}_witch_cleaned.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        witch: v\$(python3 WITCH/witch.py -v | cut -f2 -d " ")
    END_VERSIONS
    """
}
