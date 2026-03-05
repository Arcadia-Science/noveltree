process ORTHOFINDER_PHYLOHOGS {
    label 'process_medium'
    stageInMode = "copy"

    container 'arcadiascience/orthofinder_2.5.4:1.0.0'

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
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #####################################################################################
    # Prep work: set up the directory structure and file paths so OrthoFinder recognizes
    # this as a resumed run with pre-computed gene trees.

    # Get the name of the orthofinder results directory (has inflation parameter in name)
    of_results_dir=\$(ls -d Results*)
    of_working_dir=\$(grep "WorkingDirectory" \$of_results_dir/Log.txt | head -n1 | sed "s/WorkingDirectory_Base: //g")

    # Replace the original working directory paths with the current wd
    sed -i "s|\${of_working_dir}|\$(pwd)/|g" \$of_results_dir/Log.txt
    sed -i "s|OrthoFinder/Results_.*/Work|\$of_results_dir/Work|g" \$of_results_dir/Log.txt
    sed -i "s|\${of_working_dir}|\$(pwd)/|g" \$of_results_dir/WorkingDirectory/clusters_OrthoFinder_*

    # Set WorkingDirectory_Trees; OrthoFinder appends Trees_ids/ internally
    sed -i "/^WorkingDirectory_Base.*/a WorkingDirectory_Trees: \$(pwd)/\$of_results_dir/Gene_Trees/" \$of_results_dir/Log.txt

    #####################################################################################
    # Convert GeneRax trees to OrthoFinder's expected format:
    #   - Directory: Gene_Trees/Trees_ids/
    #   - Filename:  OG{7digit}_tree_id.txt
    #   - Leaf labels: OrthoFinder internal IDs (e.g. 27_153) instead of gene names
    translate_gene_trees.py SequenceIDs.txt \$of_results_dir

    #####################################################################################
    # Run orthofinder to infer hierarchical orthogroups
    orthofinder \\
        -n HOGs \\
        -s $species_tree \\
        -ft \$of_results_dir/ \\
        -a ${task.cpus} \\
        -y

    # Preserve GeneRax reconciled gene family trees in the output
    mkdir -p Results_HOGs/GeneRax_Reconciled_GFTs
    cp *_reconciled_gft.newick Results_HOGs/GeneRax_Reconciled_GFTs/

    # Clean up to avoid conflicting filenames in output
    rm -r \$of_results_dir
    mv Results_HOGs/WorkingDirectory Results_HOGs/WorkingDirectory_Hogs
    rm -f Results_HOGs/Citation.txt
    mv Results_HOGs/Log.txt Results_HOGs/Hogs_Log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$( orthofinder --help | head -n1 | sed 's/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
