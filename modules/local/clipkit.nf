process CLIPKIT {
    tag "$fasta"
    label 'process_lowcpu'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/clipkit_1.3.0-seqmagick_0.8.4:0.0.1' :
        '' }"

    publishDir(
        path: "${params.outdir}/clipkit_cleaned_msas",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(fasta)              // Filepaths to the MSAs
    val(min_ungapped_length) // Minimum ungapped length of sequences after alignment trimming

    output:
    path("*_clipkit.fa") , emit: trimmed_msas
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    # Get the name of the orthogroup we are processing
    prefix=\$(echo $fasta | cut -f1 -d "_")

    # Trim the MSAs for each orthogroup containing at least 4 species.
    clipkit ${fasta} -o \${prefix}_tmp.fa $args

    # Remove sequences with a minimum non-gapped length of 25 AA.
    seqmagick convert \\
        --min-ungapped-length $min_ungapped_length \\
        \${prefix}_tmp.fa \\
        \${prefix}_clipkit.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$( clipkit --version | sed "s/clipkit //g" )
        seqmagick: \$( seqmagick --version | cut -f2 -d" " )
    END_VERSIONS
    """
}
