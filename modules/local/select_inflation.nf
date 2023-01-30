process SELECT_INFLATION {
    tag "Selecting best MCL Inflation"
    label 'process_single'

    container "${ workflow.containerEngine == 'docker' ?
        'austinhpatton123/select_mcl_inflation_r-4.2.2_elbow_tidy_reshape_cowplot': '' }"

    publishDir(
        path: "${params.outdir}/orthogroup_summaries",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path og_summaries // Files with summaries of orthogroups inferred using a specific inflation parameter

    output:
    path "cogeqc_results.tsv" , emit: cogeqc_summary
    path "best_inflation_param.txt", emit: best_inflation
    path "inflation_summaries.pdf",  emit: summary_plot
    path "versions.yml" ,            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Pull in the summaries produced by each MCL inflation parameter and
    # identify the parameter value that produces the best quality results.
    # Combine the per-inflation-parameter tables into a single table 
    awk 'FNR==1 && NR!=1{next;}{print}' *.tsv > cogeqc_results.tsv
    
    #Clean up
    rm *_cogeqc_summary.tsv
    
    # Run the script to summarize and produce a figure of these results.
    select_inflation.R

    # And spit out the value of the selected inflation parameter to be
    # captured into a channel from stdout
    sed -i "s/\\[1] //g" best_inflation_param.txt

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        tidyverse: \$( grep "tidyverse" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        reshape: \$( grep "reshape" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        cowplot: \$( grep "cowplot" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
        elbow: \$( grep "elbow" version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}
