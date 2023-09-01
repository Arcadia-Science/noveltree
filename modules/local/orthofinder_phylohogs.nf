process ORTHOFINDER_PHYLOHOGS {
    label 'process_medium'
    stageInMode = "copy"

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder_2.5.4:1.0.0' :
        '' }"

    publishDir(
        path: "${params.outdir}/orthofinder/complete_dataset/",
        mode: params.publish_dir_mode, overwrite: false,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file species_tree        // Rooted species tree inferred using SpeciesRax
    file orthogroups         // Orthofinder results directory inferred using MCL
    file orthofinder_fastas  // Orthofinder-formatted fasta files
    file orthofinder_seq_ids // Orthofinder sequence IDs
    file orthofinder_spp_ids // Orthofinder species IDs
    file generax_gfts        // Reconciled gene family trees from GeneRax
    file blast               // Blast similarity scores

    output:
    path "Results_HOGs/" , emit: phylohogs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #####################################################################################
    # Some prep-work needs to be done prior to running orthofinder one last time,
    # parsing orthogroups into phylogenetically hierarchical orthogroups and identifying
    # orthologs. The following set of commands tidies things up so that
    # all required input files are present for orthofinder to recognize this as a
    # "resumed run", including some modification of filepaths in the orthofinder
    # log file, as these indicate where data is stored - misspecification of these
    # paths will prevent this from running.

    # Replace the original working directory with the current wd in a number of files.
    # First get the name of the orthofinder (of) results directory (will have a specific
    # inflation parameter we cannnot know a priori)
    of_results_dir=\$(ls -d Results*)
    of_working_dir=\$(grep "WorkingDirectory" \$of_results_dir/Log.txt | head -n1 | sed "s/WorkingDirectory_Base: //g")
    sed -i "s|\${of_working_dir}|\$(pwd)/|g" \$of_results_dir/Log.txt
    sed -i "s|OrthoFinder/Results_.*/Work|\$of_results_dir/Work|g" \$of_results_dir/Log.txt
    sed -i "s|\${of_working_dir}|\$(pwd)/|g" \$of_results_dir/WorkingDirectory/clusters_OrthoFinder_*

    # Now, add the "WorkingDirectory_Trees" path to trick orthofinder into recognizing the input
    sed -i "/^WorkingDirectory_Base.*/a WorkingDirectory_Trees: \$(pwd)/\$of_results_dir/Gene_Trees/" \$of_results_dir/Log.txt

    # And move the directory containing gene trees (currently in the working directory) within the results dir
    mkdir \$of_results_dir/Gene_Trees && mv *gft.newick \$of_results_dir/Gene_Trees/

    #####################################################################################

    # Run orthofinder to sort into hierarchical orthogoups.
    orthofinder \
    -n HOGs \
    -s $species_tree \
    -ft \$of_results_dir/ \
    -a ${task.cpus} \
    -y

    # Move the generax reconciled gene family trees within the the Orthofinder HOG directory
    mv \$of_results_dir/Gene_Trees/ Results_HOGs/GeneRax_Reconciled_GFTs

    # And clean up,rename a few things so as not to have conflicting filenames in the resultant output
    rm -r \$of_results_dir
    mv Results_HOGs/WorkingDirectory Results_HOGs/WorkingDirectory_Hogs
    rm Results_HOGs/Citation.txt
    mv Results_HOGs/Log.txt Results_HOGs/Hogs_Log.txt
    """
}
