process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.3--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.3--hdfd78af_0' }"

    input:
    each mcl_inflation
    val ch_similarities

    output:
    path "MCL-*-fpath.txt" , emit: og_fpath

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inflation_param = mcl_inflation ? "${mcl_inflation}" : '1.5'
    """
    # add a sleep time equal to the inflation parameter in case things get overloaded when running locally?
    sleep $inflation_param
    orthofinder \\
        -b ../../../${params.outdir}/orthofinder/WorkingDirectory/ \\
        -n "Inflation_${inflation_param}" \\
        -I $inflation_param \\
        -M msa -X -os -z \\
        -a ${task.cpus} \\
        $args
    
    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mv ../../../${params.outdir}/orthofinder/WorkingDirectory/OrthoFinder/Results_Inflation_${inflation_param}/ ../../../${params.outdir}/orthofinder/
    
    # And output a file used to track job completion and specify filepaths downstream
    ( cd ../../../${params.outdir}/orthofinder/Results_Inflation_${inflation_param}/; pwd ) > MCL-${inflation_param}-fpath.txt
    
    """
}

