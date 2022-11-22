process CLIPKIT {
    tag "MSA trimming"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bioconductor-uniprot.ws==2.34.0--r41hdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-uniprot.ws:2.34.0--r41hdfd78af_0':
        'docker.io/austinhpatton123/clipkit' }"
        
    publishDir(
        path: "${params.outdir}/trimmed-msas",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)   // Filepaths to the MSAs with at least four species 

    output:
    tuple val(meta), path("*-clipkit.fa")      , emit: trimmed_msas
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    # Trim the MSAs for each orthogroup containing at least 4 species. 
    clipkit ${fasta} -o ${meta.og}-clipkit.fa
    
    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    clipkit: \$( cat version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    #END_VERSIONS
    """
}

