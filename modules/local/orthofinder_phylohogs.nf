process ORTHOFINDER_PHYLOHOGS {
    // label 'process_highthread'
    label 'process_medium'
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/orthofinder-2.5.4_r-4.2.2' : 
        '' }"

    publishDir(
        path: "${params.outdir}/orthofinder_phylohogs/",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
    input:
    path species_tree    // Rooted species tree inferred using SpeciesRax
    path speciesrax_gfts // Reconciled gene family trees used in species tree inference
    path generax_gfts    // Remaining reconciled gene family trees
    

    output:
    path "*" , emit: phylohogs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Copy the previous orthofinder results into the current directory
    cp -r ${params.outdir}/orthofinder/complete_dataset/* .
    
    # Rename it
    mv Results_Inflation_*/ orthogroups
    
    # Replace the original working directory with the current wd in a number of files. 
    dir=\$(grep "WorkingDirectory" orthogroups/Log.txt | head -n1 | sed "s/WorkingDirectory_Base: //g")
    sed -i "s|\${dir}|\$(pwd)/|g" orthogroups/Log.txt
    sed -i "s|OrthoFinder/Results_.*/Work|orthogroups/Work|g" orthogroups/Log.txt
    sed -i "s|\${dir}|\$(pwd)/|g" orthogroups/WorkingDirectory/clusters_OrthoFinder_*
    sed -i "s/-test-proteome.fasta//g" orthogroups/Log.txt # Delete eventually - proteomes should be named exactly what we want tip labels to be. 
    
    # Now, add the "WorkingDirectory_Trees" path to trick orthofinder into recognizing the input 
    sed -i "/^WorkingDirectory_Base.*/a WorkingDirectory_Trees: \$(pwd)/orthogroups/Gene_Trees/" orthogroups/Log.txt
    
    # Now make the directory, and populate with reconciled gene trees
    mkdir orthogroups/Gene_Trees
    
    # And copy all of the gene family trees into that directory. 
    mv OG*.newick orthogroups/Gene_Trees
    
    # And bring the orthofinder formatted Blast data in 
    cp ${params.outdir}/diamond/Blast* .
    
    # Run orthofinder to sort into hierarchical orthogoups. 
    orthofinder \
    -n HOGs \
    -s $species_tree \
    -ft orthogroups/ \
    -a 60 -y 
    
    # Lastly, clean up all the duplicated files that are in the current directory, 
    # and then move all orthofinder output to the current working directory 
    # to clean up
    rm Blast*
    rm diamond*
    rm *.newick
    rm *.txt
    rm *.fa
    rm -r orthogroups

    mv Results_HOGs/* . && rm -r Results_HOGs/
    """
}

