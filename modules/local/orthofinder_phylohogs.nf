process ORTHOFINDER_PHYLOHOGS {
    // label 'process_highthread'
    label 'process_medium'
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/orthofinder-2.5.4_r-4.2.2' : 
        '' }"
    stageInMode = "copy"
    publishDir(
        path: "${params.outdir}/orthofinder_phylohogs/",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
    input:
    path species_tree                              // Rooted species tree inferred using SpeciesRax
    path orthogroups   // Orthofinder results directory inferred using MCL
    path orthofinder_fastas //
    path orthofinder_seq_ids //
    path orthofinder_spp_ids //
    path speciesrax_gfts   // Reconciled gene family trees used in species tree inference
    path generax_gfts    // Remaining reconciled gene family trees
    path blast                                     // Blast similarity scores

    output:
    path "*" , emit: phylohogs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Replace the original working directory with the current wd in a number of files. 
    # First get the name of the Results directory (will have a specific inflation parameter we cannnot know a priori)
    results_dir=\$(ls -d Results*)
    dir=\$(grep "WorkingDirectory" \$results_dir/Log.txt | head -n1 | sed "s/WorkingDirectory_Base: //g")
    sed -i "s|\${dir}|\$(pwd)/|g" \$results_dir/Log.txt
    sed -i "s|OrthoFinder/Results_.*/Work|\$results_dir/Work|g" \$results_dir/Log.txt
    sed -i "s|\${dir}|\$(pwd)/|g" \$results_dir/WorkingDirectory/clusters_OrthoFinder_*
    sed -i "s/-test-proteome.fasta//g" \$results_dir/Log.txt # Delete eventually - proteomes should be named exactly what we want tip labels to be. 
    
    # Now, add the "WorkingDirectory_Trees" path to trick orthofinder into recognizing the input 
    sed -i "/^WorkingDirectory_Base.*/a WorkingDirectory_Trees: \$(pwd)/\$results_dir/Gene_Trees/" \$results_dir/Log.txt
    
    # And move the directory containing gene trees (currently in the working directory) within the results dir
    mkdir \$results_dir/Gene_Trees && mv *gft.newick \$results_dir/Gene_Trees/

    # Now make the directory, and populate with reconciled gene trees
    # mkdir Results_Inflation_*/Gene_Trees
    
    # And copy all of the gene family trees into that directory. 
    # mv OG*.newick orthogroups/Gene_Trees
    
    # Run orthofinder to sort into hierarchical orthogoups. 
    orthofinder \
    -n HOGs \
    -s $species_tree \
    -ft \$results_dir/ \
    -a ${task.cpus} \
    -y
    """
}

