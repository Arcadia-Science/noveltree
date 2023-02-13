process SPECIESRAX {
    tag "SpeciesRax"
    //label 'process_highthread' // Possible specification for full analysis
    label 'process_medium' // Used for debugging

    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/generax:19604b7': '' }"

    publishDir(
        path: "${params.outdir}/speciesrax",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path init_species_tree // Filepath to the starting species tree
    path generax_map       // Filepath to the generax gene-species map file
    path gene_trees        // Filepaths to the starting gene trees
    path alignments        // Filepaths to the gene family alignments
    path families          // Filepath to the families file

    output:
    path "*"                                          , emit: results
    path "species_trees/inferred_species_tree.newick" , emit: speciesrax_tree
    path "**_reconciled_gft.newick"                   , emit: speciesrax_gfts
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $init_species_tree \\
        --families $families \\
        --prefix SpeciesRax \\
        $args

    # Remove the redundant result directory, moving everything into the
    # working directory
    mv SpeciesRax/* .
    rm -r SpeciesRax

    # Rename the inferred reconciled gene trees to be named after their corresponding orthogroup
    for og in \$(ls results/)
    do
        mv results/\$og/*.newick results/\$og/\${og}_reconciled_gft.newick
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$( generax | head -n1 | sed "s/.*GeneRax //g" )
    END_VERSIONS
    """
}
