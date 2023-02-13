process ORTHOFINDER_PHYLOHOGS {
    // label 'process_highthread'
    label 'process_medium'
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder:2.5.4' :
        '' }"
    stageInMode = "copy"
    publishDir(
        path: "${params.outdir}/orthofinder_phylohogs/",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path species_tree        // Rooted species tree inferred using SpeciesRax
    path orthogroups         // Orthofinder results directory inferred using MCL
    path orthofinder_fastas  // Orthofinder-formatted fasta files
    path orthofinder_seq_ids // Orthofinder sequence IDs
    path orthofinder_spp_ids // Orthofinder species IDs
    path speciesrax_gfts     // Reconciled gene family trees used in species tree inference
    path generax_gfts        // Remaining reconciled gene family trees
    path blast               // Blast similarity scores

    output:
    path "*" , emit: phylohogs

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
    """
}
