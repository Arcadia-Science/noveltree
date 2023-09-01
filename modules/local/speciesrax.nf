process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_generax'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems.

    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/generax_19604b71:1.0.0': '' }"

    publishDir(
        path: "${params.outdir}/speciesrax",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file map_links       // Filepath to the generax gene-species map file
    file gene_trees      // Filepaths to the starting gene trees
    file alignments      // Filepaths to the gene family alignments
    file rooted_spp_tree // Filepath to the rooted asteroid species tree

    output:
    path "*"                                          , emit: results
    path "species_trees/inferred_species_tree.newick" , emit: speciesrax_tree
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def starting_tree = (rooted_spp_tree && file(rooted_spp_tree).exists()) ? rooted_spp_tree : 'MiniNJ'
    """
    # Recode selenocysteine as a gap character:
    # RAxML-NG (used under the hood by SpeciesRax and
    # GeneRax) cannot handle these. Even if rare,
    # their inclusion leads a number of gene families
    # to be excluded from analyses.
    sed -E -i '/>/!s/U/-/g' *.fa

    # Do the same for Pyrrolysine
    sed -E -i '/>/!s/O/-/g' *.fa

    # Construct the family files for each gene family
    echo "[FAMILIES]" > speciesrax_orthogroup.families
    for msa in \$(ls *fa)
    do
        # Get the OG name
        og=\$(echo \$msa | cut -f1 -d"_")
        tree=\$(ls \${og}*.treefile)

        # Populate the families file for this gene family for the
        # analysis with SpeciesRax
        # We will be using LG+G4+F for all gene families
        echo "- \${og}" >> speciesrax_orthogroup.families
        echo "starting_gene_tree = \${tree}" >> speciesrax_orthogroup.families
        echo "mapping = \${og}_map.link" >> speciesrax_orthogroup.families
        echo "alignment = \$msa" >> speciesrax_orthogroup.families
        echo "subst_model = LG+G4+F" >> speciesrax_orthogroup.families
    done


    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $starting_tree \\
        --families speciesrax_orthogroup.families \\
        --prefix SpeciesRax \\
        --strategy SKIP \\
        --si-estimate-bl \\
        --per-species-rates \\
        $args

    # Remove the redundant result directory, moving everything into the
    # working directory, deleiting the meaningless reconciliations
    # directory and cleaning up
    mv SpeciesRax/* .
    rm -r reconciliations
    rm -r SpeciesRax

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$( generax | head -n1 | sed "s/.*GeneRax //g" )
    END_VERSIONS
    """
}
