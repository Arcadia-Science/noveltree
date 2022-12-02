process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_highthread'
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/orthofinder-2.5.4_r-4.2.2' : 
        '' }"

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
        -a ${task.cpus} \\
        $args

    # Clean up the input files from the working directory so they're not returned. 
    rm ./*.dmnd 
    rm ./*.fa 
    rm ./*IDs.txt
    
    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mv ./OrthoFinder/Results_Inflation_${inflation_param}/ \$dataDir/

    # And output a file used to track job completion and specify filepaths downstream
    ( cd \$dataDir/Results_Inflation_${inflation_param}/; pwd ) > MCL-${inflation_param}-fpath.txt

    """
}

