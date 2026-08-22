process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_generax'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems.

    container 'arcadiascience/generax_56f3ed0:1.1.3'

    storeDir "${params.outdir}/species_trees/speciesrax"

    input:
    file input_bundles   // Sharded archives of trees, maps, and manifests
    file rooted_spp_tree // Filepath to the rooted asteroid species tree
    val species_names    // Complete list of input-dataset species names

    output:
    path "inferred_species_tree.newick"  , emit: speciesrax_tree
    path "starting_species_tree.newick"
    path "species_tree_*.newick"
    path "*.txt"
    path "generax.log"
    path "speciesrax_orthogroup.families"
    path "speciesrax_gene_tree_validation.tsv"
    path "speciesrax_selected_families.tsv"
    path "speciesrax_family_selection.tsv"
    path "speciesrax_selected_species_coverage.tsv"

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

    # SpeciesRax requires rooted-binary Newick input, while an otherwise binary
    # unrooted tree is commonly serialized with a three-child root. Resolve all
    # multifurcations deterministically in these task-local copies, then fail
    # before MPI initialization if any tree is still not strictly binary.
    # Stored upstream gene trees are never modified.
    prepare_speciesrax_gene_trees.py \
        --manifest-glob 'speciesrax_inputs_*.tsv' \
        --report speciesrax_gene_tree_validation.tsv

    # Preserve the complete run-level denominator for occupancy and leaf caps.
    echo "$species_names" \
        | sed "s/\\[//g; s/\\]//g" \
        | tr "," "\\n" \
        | sed "s/^[[:space:]]*//; s/[[:space:]]*\$//" \
        > speciesrax_expected_species.txt

    # Temporarily enforce the production-scale SpeciesRax limits here, after
    # final tree validation. These cutoffs intentionally retain multicopy
    # families. Once validated across datasets, move the same classification
    # into upstream orthogroup routing to avoid staging rejected families.
    select_speciesrax_families.py \
        --manifest-glob 'speciesrax_inputs_*.tsv' \
        --validation-report speciesrax_gene_tree_validation.tsv \
        --expected-species-file speciesrax_expected_species.txt \
        --min-species-occupancy ${params.speciesrax_min_species_occupancy} \
        --max-mean-copies ${params.speciesrax_max_mean_copies} \
        --max-copies-per-species ${params.speciesrax_max_copies_per_species} \
        --max-total-leaves-factor ${params.speciesrax_max_total_leaves_factor} \
        --selected-manifest speciesrax_selected_families.tsv \
        --report speciesrax_family_selection.tsv \
        --species-coverage-report speciesrax_selected_species_coverage.tsv

    # Construct the family file from the selected, validated families.
    echo "[FAMILIES]" > speciesrax_orthogroup.families
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
    done < speciesrax_selected_families.tsv
    rm -f speciesrax_inputs_*.tsv speciesrax_expected_species.txt


    # Do not request per-species or per-family rates: the REROOT path uses the
    # shared global D/L parameterization.
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
