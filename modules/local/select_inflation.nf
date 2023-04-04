process SELECT_INFLATION {
    tag "Selecting best MCL Inflation"
    label 'process_low'

    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/select_mcl_inflation_params:0.0.1': '' }"

    publishDir(
        path: "${params.outdir}/orthogroup_summaries",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file og_summaries // Files with summaries of orthogroups inferred using a specific inflation parameter
    val min_spp       // Minimum number of species for orthogroup phylogenetic inference.

    output:
    path "cogeqc_results.tsv"       , emit: cogeqc_summary
    path "best_inflation_param.txt" , emit: best_inflation
    path "inflation_summaries.pdf"  , emit: summary_plot
    path "versions.yml"             , emit: versions

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
    select_inflation.R ${min_spp}
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
