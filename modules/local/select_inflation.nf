process SELECT_INFLATION {
    tag "Selecting best MCL Inflation"
    label 'process_single'

    container "${ workflow.containerEngine == 'docker' ?
        'austinhpatton123/select_mcl_inflation_r-4.2.2_elbow_tidy_reshape_cowplot': '' }"
        
    publishDir(
        path: "${params.outdir}/orthogroup-summaries",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path og_summaries // Files with summaries of orthogroups inferred using a specific inflation parameter

    output:
    path "best-inflation-param.txt", emit: best_inflation
    path "inflation_summaries.pdf",  emit: summary_plot
    path "versions.yml" ,            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Pull in the summaries produced by each MCL inflation parameter and 
    # identify the parameter value that produces the best quality results.

    # Run the script to summarize and produce a figure of these results. 
    Rscript $projectDir/bin/select_inflation.R
    
    # And spit out the value of the selected inflation parameter to be 
    # captured into a channel from stdout 
    sed -i "s/\\[1] //g" best-inflation-param.txt
    
    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: \$( grep "tidyverse" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        reshape: \$( grep "reshape" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        cowplot: \$( grep "cowplot" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        elbow: \$( grep "elbow" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}