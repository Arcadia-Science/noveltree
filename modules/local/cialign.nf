process CIALIGN {
    tag "$fasta"
    label 'process_lowcpu'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/cialign_1.1.0:0.0.1' :
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/cialign_1.1.0:0.0.1' :
        '' }"

    publishDir(
        path: "${params.outdir}/trimmed_msas",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(fasta)   // Filepaths to the MSAs

    output:
    path("*_cialign.fa") , emit: trimmed_msas
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    # Get the name of the orthogroup we are processing
    prefix=\$(echo ${fasta} | sed "s/_.*//g")

    # Clean up the MSAs for each orthogroup containing at least 4 species.
    python CIAlign_1.1.0/CIAlign/CIAlign.py \
        --infile ${fasta} \
        --outfile_stem="\${prefix}" \
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CIAlign: \$( CIAlign --version )
    END_VERSIONS
    """
}
