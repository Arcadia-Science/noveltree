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
    path "MCL-*-fpath.txt"         , emit: og_fpath
    path "four_spp_ogs.txt"        , emit: four_spp_ogs

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
    #mv ./OrthoFinder/Results_Inflation_${inflation_param}/ \$dataDir/
    #mv ../../../${params.outdir}/orthofinder/WorkingDirectory/OrthoFinder/Results_Inflation_${inflation_param}/ ../../../${params.outdir}/orthofinder/
    
    # And output a file used to track job completion and specify filepaths downstream
    #( cd ../../../${params.outdir}/orthofinder/Results_Inflation_${inflation_param}/; pwd ) > MCL-${inflation_param}-fpath.txt
    ( cd \$dataDir/Results_Inflation_${inflation_param}/; pwd ) > MCL-${inflation_param}-fpath.txt

    # Identify the four species OGs if this is the full analysis
    ogSppCounts=\$dataDir/Results_Inflation_${inflation_param}/Orthogroups/Orthogroups.GeneCount.tsv
    
    # Convert the gene counts to a binary
    tail -n+2 \$ogSppCounts | awk '{for (i=2; i<= NF; i++) {if(\$i >= 1) { \$i=1; }}print }' > tmp
    
    # Pull out the orthogroup names and the number of species there's in
    cut -f1 -d" " \$ogSppCounts > og_names.txt
    cut -f2- tmp | awk 'NF{NF-=1};1' | awk '{sum=0; for(i=1; i<=NF; i++) sum += \$i; print sum}' > counts.txt
    
    # Combine them
    paste og_names.txt counts.txt > og_spp_counts.txt
    awk '(NR>1) && (\$2 >= 9 ) ' og_spp_counts.txt | cut -f1 > four_spp_ogs.txt
    
    # Now, get the corresponding filepaths
    msaDir=\$( cd \$dataDir/Results_Inflation_${inflation_param}/Orthogroup_Sequences/; pwd )
    sed "s|.*|\${msaDir}/&.fa|g" four_spp_ogs.txt > four_spp_og_fpaths.txt
    
    # combine these into what will be turned into a groovy map and file input. 
    paste -d"," four_spp_ogs.txt four_spp_og_fpaths.txt > tmp
    echo "orthogroup,file" > header
    cat header tmp > four_spp_ogs.txt
    """
}

