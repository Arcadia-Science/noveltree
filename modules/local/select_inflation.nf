process SELECT_INFLATION {
    tag "Selecting best MCL Inflation"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bioconductor-r-tidyverse==1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1':
        'austinhpatton123/select_inflation_param:latest' }"
        
    publishDir(
        path: "${params.outdir}/orthogroup-summaries",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path og_summaries // Files with summaries of orthogroups inferred using a specific inflation parameter

    output:
    path "best-inflation-param.txt", emit: best_inflation
    path "inflation_summaries.pdf", emit: summary_plot
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Pull in the summaries produced by each MCL inflation parameter and 
    # identify the parameter value that produces the best quality results.
        
    # Start by assembling all summaries into a single table
    #cat *summary.tsv | head -n1 > inflation_summaries.tsv # Pull out the header
    
    # Now add the summaries
    #cat *summary.tsv | grep -v "Inflation" >> inflation_summaries.tsv # Pull out the header

    # Run the script to summarize and produce a figure of these results. 
    Rscript $projectDir/bin/select_inflation.R
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: \$( grep "tidyverse" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        reshape: \$( grep "reshape" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        cowplot: \$( grep "cowplot" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        elbow: \$( grep "elbow" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}