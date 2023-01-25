process CLIPKIT {
    tag "$fasta"
    label 'process_medium'

    container "${ workflow.containerEngine == 'docker' ? 'docker.io/austinhpatton123/clipkit' :
        '' }"

    publishDir(
        path: "${params.outdir}/trimmed_msas",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(fasta)   // Filepaths to the MSAs

    output:
    path("*_clipkit.fa")                  , emit: trimmed_msas
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    prefix=\$(basename "${fasta}" _mafft.fa)

    # Trim the MSAs for each orthogroup containing at least 4 species.
    clipkit ${fasta} -o \${prefix}_clipkit.fa $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$( clipkit --version | sed "s/clipkit //g" )
    END_VERSIONS
    """
}
