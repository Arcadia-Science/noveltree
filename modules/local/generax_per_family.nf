process GENERAX_PER_FAMILY {
    tag "GeneRax: Per-family rates"
    label 'process_generax_per_family'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems. 
    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/generax_19604b7:0.0.1': '' }"

    publishDir(
        path: "${params.outdir}/generax/per_family_rates",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    // The following tuple includes (for each gene family) the 
    // gene-species map files, gene trees, family files, and the 
    // species tree inferred from SpeciesRax
    tuple file(generax_map), file(gene_tree), file(alignment), file(family), file(species_tree)

    output:
    path "*"                        , emit: results
    path "**_reconciled_gft.newick" , emit: generax_per_fam_gfts

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

    # Get the name of the focal gene family
    og=\$(echo ${family} | sed "s/.family//g")

    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $species_tree \\
        --families $family \\
        --per-family-rates \\
        --prefix \$og \\
        $args

    # Clean up
    rm -r \$og/gene_optimization_*

    # Rename the inferred reconciled gene trees to be named after their corresponding orthogroup
    mv \$og/results/\$og/geneTree.newick \$og/results/\$og/\${og}_reconciled_gft.newick
    """
}
