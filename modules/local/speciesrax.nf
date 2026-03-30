process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_generax'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems.

    container 'arcadiascience/generax_56f3ed0:1.1.3'

    storeDir "${params.outdir}/species_trees/speciesrax"

    input:
    file map_links       // Filepath to the generax gene-species map file
    file gene_trees      // Filepaths to the starting gene trees
    file alignments      // Filepaths to the gene family alignments
    file rooted_spp_tree // Filepath to the rooted asteroid species tree

    output:
    path "inferred_species_tree.newick"  , emit: speciesrax_tree
    path "starting_species_tree.newick"
    path "species_tree_*.newick"
    path "*.txt"
    path "generax.log"
    path "speciesrax_orthogroup.families"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def starting_tree = (rooted_spp_tree && file(rooted_spp_tree).exists()) ? rooted_spp_tree : "MiniNJ"
    """
    # Construct the family files for each gene family
    echo "[FAMILIES]" > speciesrax_orthogroup.families
    for msa in \$(ls *fa)
    do
        # Get the OG name
        og=\$(echo \$msa | cut -f1 -d"_")
        tree=\$(ls \${og}*.newick)
        map_link=\$(ls \${og}*_map.link | head -n1)

        # Populate the families file for this gene family for the
        # analysis with SpeciesRax
        # We will be using LG+G4+F for all gene families
        echo "- \${og}" >> speciesrax_orthogroup.families
        echo "starting_gene_tree = \${tree}" >> speciesrax_orthogroup.families
        echo "mapping = \${map_link}" >> speciesrax_orthogroup.families
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

    # Move SpeciesRax output into the working directory and clean up
    mv SpeciesRax/* .
    rm -rf reconciliations results SpeciesRax

    # Flatten species_trees/ into working directory
    mv species_trees/* .
    rm -r species_trees

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$( generax | head -n1 | sed "s/.*GeneRax //g" )
    END_VERSIONS
    """
}
