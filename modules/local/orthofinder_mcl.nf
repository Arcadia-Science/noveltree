process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_mcl'
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/orthofinder-2.5.4_r-4.2.2' : 
        '' }"

    input:
    each mcl_inflation
    val ch_similarities
    val mcl_test
    val min_num_spp
    val min_num_groups
    val max_copy_num_filt1
    val max_copy_num_filt2

    output:
    path "MCL-*-fpath.txt"               , emit: og_fpath
    path "all_ogs_counts.csv"            , emit: all_ogs
    path "extreme_core_ogs_counts.csv"   , emit: extreme_core_ogs
    path "remaining_core_ogs_counts.csv" , emit: remaining_core_ogs
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inflation_param = mcl_inflation ? "${mcl_inflation}" : '1.5'
    def num_spp = min_num_spp ? "${min_num_spp}" : '4'
    def num_grp = min_num_groups ? "${min_num_groups}" : '2'
    def copy_num1 = max_copy_num_filt1 ? "${max_copy_num_filt1}" : '5'
    def copy_num2 = max_copy_num_filt2 ? "${max_copy_num_filt2}" : '10'
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

    # The below does some reformatting to prepare for MSA inference
    # in the case that we're doing a full analysis.
    # If instead we're just doing the inflation parameter testing,
    # output some filler files to save compute
    if [ "$testing_mcl" == "true" ]; then
        echo "Testing mcl inflation parameters." > all_ogs_counts.csv
        echo "Testing mcl inflation parameters." > extreme_core_ogs_counts.csv
        echo "Testing mcl inflation parameters." > remaining_core_ogs_counts.csv
    else
        # Running a full analysis - go ahead and pull out these OGs.
        # Identify the core OGs if this is the full analysis
        ogSppCounts=\$dataDir/Results_Inflation_${inflation_param}/Orthogroups/Orthogroups.GeneCount.tsv
        
        # Run the scripts to generate the orthogroup species/taxa gene count summaries and filtered sets
        Rscript $projectDir/bin/og-tax-summary.R \$ogSppCounts ${params.input} $num_spp $num_grp $copy_num1 $copy_num2
        
        # Add a column of filepaths to these
        msaDir=\$( cd \$dataDir/Results_Inflation_${inflation_param}/Orthogroup_Sequences/; pwd )
        # Create the column
        echo "file" > colname.txt
        tail -n+2 all_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msaDir}/&.fa|g" | sed '1 i\\file' > all_og_fpaths.txt
        tail -n+2 extreme_core_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msaDir}/&.fa|g" | sed '1 i\\file' > extr_core_og_fpaths.txt
        tail -n+2 remaining_core_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msaDir}/&.fa|g" | sed '1 i\\file' > rem_core_og_fpaths.txt
        
        # Combine them
        paste -d"," all_ogs_counts.csv all_og_fpaths.txt > tmp && mv tmp all_og_fpaths.csv
        paste -d"," extreme_core_ogs_counts.csv extr_core_og_fpaths.txt > tmp && mv tmp extreme_core_ogs_counts.csv
        paste -d"," remaining_core_ogs_counts.csv rem_core_og_fpaths.txt > tmp && mv tmp remaining_core_ogs_counts.csv
        
        # Clean up
        rm *og_fpaths.txt
    fi   
    """
}

