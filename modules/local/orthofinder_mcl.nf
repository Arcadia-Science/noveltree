process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_mcl'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.3--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.3--hdfd78af_0' }"

    input:
    each mcl_inflation
    val ch_similarities
    val mcl_test
    val min_num_spp
    val max_copy_num

    output:
    path "MCL-*-fpath.txt"      , emit: og_fpath
    path "core_ogs*.csv"        , emit: core_ogs
    path "all_ogs*.csv"         , emit: all_ogs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inflation_param = mcl_inflation ? "${mcl_inflation}" : '1.5'
    def num_spp = min_num_spp ? "${min_num_spp}" : '2'
    def copy_num = max_copy_num ? "${max_copy_num}" : '50'
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

    # The below does some reformatting to prepare for MSA inference
    # in the case that we're doing a full analysis.
    # If instead we're just doing the inflation parameter testing,
    # output some filler files to save compute
    if [ "$testing_mcl" == "true" ]; then
        echo "Testing mcl inflation parameters." > all_ogs_test_${inflation_param}.csv
        echo "Testing mcl inflation parameters." > core_ogs_test_${inflation_param}.csv
    else
        # Running a full analysis - go ahead and pull out these OGs.
        # Identify the core OGs if this is the full analysis
        ogSppCounts=\$dataDir/Results_Inflation_${inflation_param}/Orthogroups/Orthogroups.GeneCount.tsv
        
        # Convert the gene counts to a binary
        tail -n+2 \$ogSppCounts | awk '{for (i=2; i<= NF; i++) {if(\$i >= 1) { \$i=1; }}print }' > tmp
        
        # Pull out the orthogroup names and the number of species there's in
        cut -f1 -d" " tmp > og_names.txt
        cut -f2- -d" " tmp | awk 'NF{NF-=1};1' | awk '{sum=0; for(i=1; i<=NF; i++) sum += \$i; print sum}' > counts.txt
    
        # Combine them
        paste -d"," og_names.txt counts.txt > og_spp_counts.txt
        
        # Get the total number of gene copies
        tail -n+2 \$ogSppCounts | awk '{print \$NF}' > total_counts.txt
        
        # and from this, the per-species mean copy number (for spp with the OG).
        paste total_counts.txt counts.txt | awk '{print(\$1/\$2)}' > spp_mean_counts.txt
        
        # combine with other metadata
        paste -d"," og_spp_counts.txt spp_mean_counts.txt > tmp && mv tmp og_spp_counts.txt
    
        # Generate the filepaths to each
        msaDir=\$( cd \$dataDir/Results_Inflation_${inflation_param}/Orthogroup_Sequences/; pwd )
        sed "s|.*|\${msaDir}/&.fa|g" og_names.txt > og_fpaths.txt
    
        # combine these into what will be turned into a groovy map and file input. 
        paste -d"," og_spp_counts.txt og_fpaths.txt > tmp
        echo "orthogroup,num_spp,mean_copy_num,file" > header
        cat header tmp > all_ogs_final_${inflation_param}.csv
        
        # Pull out the OGs with at least the specified number of species with at most the specified mean copy number per species
        awk -F"," '(NR>=1) && (\$2 >= $num_spp ) && (\$3 <= $copy_num )' all_ogs_final_${inflation_param}.csv > tmp
        cat header tmp > core_ogs_final_${inflation_param}.csv
    fi   
    """
}

