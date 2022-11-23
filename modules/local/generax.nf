process GENERAX {
    tag "meta.og"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::generax==2.0.4--h19e7193_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/generax:2.0.4--h19e7193_0':
        'quay.io/biocontainers/generax:2.0.4--h19e7193_0' }"
        
    publishDir(
        path: "${params.outdir}/trimmed-msas",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)   // Filepaths to the MSAs with at least four species 

    output:
    tuple val(meta), path("*-clipkit.fa") , emit: trimmed_msas
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    # Trim the MSAs for each orthogroup containing at least 4 species. 
    clipkit ${fasta} -o ${meta.og}-clipkit.fa
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$( clipkit --version | sed "s/clipkit //g" )
    END_VERSIONS
    """
}

