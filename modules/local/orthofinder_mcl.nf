process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.3--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.3--hdfd78af_0' }"

    input:
    each mcl_inflation
    val ch_similarities
    val mcl_test

    output:
    path "MCL-*-fpath.txt" , emit: og_fpath

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inflation_param = mcl_inflation ? "${mcl_inflation}" : '1.5'
    def testing_mcl = mcl_test.equals('true') ? "${mcl_test}" : "false"
    """
    # add a sleep time equal to the inflation parameter in case things get overloaded when running locally?
    sleep \$(( \$RANDOM % 10 + 1 ))

    # The following depends on whether we're testing mcl or not. 
    if [ "$testing_mcl" == "true" ]; then
        dataDir="../../../${params.outdir}/orthofinder/mcl_testing"
        cp ../../../${params.outdir}/diamond/TestBlast* .
        for f in \$(ls TestBlast*)
        do
            mv \$f \$(echo \$f | sed "s/TestBlast/Blast/g")
        done
    else
        dataDir="../../../${params.outdir}/orthofinder/full_analysis"
        cp ../../../${params.outdir}/diamond/Blast* .
    fi
    
    cp \$dataDir/data/* .

    orthofinder \\
        -b ./ \\
        -n "Inflation_${inflation_param}" \\
        -I $inflation_param \\
        -M msa -X -os -z \\
        -a 30 \\
        $args
    # -a ${task.cpus} \\

    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mv ./OrthoFinder/Results_Inflation_${inflation_param}/ \$dataDir/
    #mv ../../../${params.outdir}/orthofinder/WorkingDirectory/OrthoFinder/Results_Inflation_${inflation_param}/ ../../../${params.outdir}/orthofinder/
    
    # And output a file used to track job completion and specify filepaths downstream
    #( cd ../../../${params.outdir}/orthofinder/Results_Inflation_${inflation_param}/; pwd ) > MCL-${inflation_param}-fpath.txt
    ( cd \$dataDir/Results_Inflation_${inflation_param}/; pwd ) > MCL-${inflation_param}-fpath.txt

    # If this is the full analysis, do some additional cleanup
    #if [ "$testing_mcl" == "false" ]; then
    #    rmdir ../../../${params.outdir}/orthofinder/WorkingDirectory/OrthoFinder/ || true
    #fi

    """
}

