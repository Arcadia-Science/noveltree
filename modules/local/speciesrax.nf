process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_generax'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems.

    container 'arcadiascience/generax_56f3ed0:1.1.3'

    storeDir "${params.outdir}/species_trees/speciesrax"

    input:
    file input_bundles   // Sharded archives of trees, maps, and manifests
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
    # Extract uncompressed shards locally. The archive layer reduces S3/API
    # fan-in without spending CPU recompressing files GeneRax immediately reads.
    for archive in speciesrax_inputs_*.tar
    do
        tar -xf "\$archive"
        rm -f "\$archive"
    done

    # Construct the family file from the validated shard manifests.
    echo "[FAMILIES]" > speciesrax_orthogroup.families
    for manifest in speciesrax_inputs_*.tsv
    do
        while IFS=\$'\t' read -r og tree map_link
        do
            if [ "\$og" = "orthogroup" ]; then
                continue
            fi
            for required in "\$tree" "\$map_link"
            do
                if [ ! -s "\$required" ]; then
                    echo "Missing or empty SpeciesRax input: \$required" >&2
                    exit 1
                fi
            done

            echo "- \${og}" >> speciesrax_orthogroup.families
            echo "starting_gene_tree = \${tree}" >> speciesrax_orthogroup.families
            echo "mapping = \${map_link}" >> speciesrax_orthogroup.families
        done < "\$manifest"
    done
    rm -f speciesrax_inputs_*.tsv


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
