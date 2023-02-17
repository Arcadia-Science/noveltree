process GENERAX {
    tag "GeneRax"
    //label 'process_highthread' // Possible specification for full analysis
    label 'process_medium' // Used for debugging

    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/generax:19604b7': '' }"

    publishDir(
        path: "${params.outdir}/generax",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path init_species_tree // Filepath to the SpeciesRax species tree
    path generax_map       // Filepath to the generax gene-species map file
    path gene_trees        // Filepaths to the starting gene trees
    path alignments        // Filepaths to the gene family alignments
    path families          // Filepath to the families file

    output:
    path "*"                        , emit: results
    path "**_reconciled_gft.newick" , emit: generax_gfts

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''

    """
    # Recode selenocysteine as a gap character:
    # RAxML-NG (used under the hood by SpeciesRax and 
    # GeneRax) cannot handle these. Even if rare, 
    # their inclusion leads a number of gene families 
    # to be excluded from analyses. 
    sed -E -i '/>/!s/U/-/g' *.fa 
    
    # Do the same for Pyrrolysine
    sed -E -i '/>/!s/O/-/g' *.fa 

    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $init_species_tree \\
        --families $families \\
        --prefix GeneRax \\
        $args

    # Remove the redundant result directory, moving everything into the
    # working directory
    mv GeneRax/* .
    rm -r GeneRax

    # Rename the inferred reconciled gene trees to be named after their corresponding orthogroup
    for og in \$(ls results/)
    do
        mv results/\$og/*.newick results/\$og/\${og}_reconciled_gft.newick
    done
    """
}
