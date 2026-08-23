process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_generax'
    stageInMode 'copy'

    // 1.1.4 adds Biopython for deterministic in-process root transfer.
    container 'arcadiascience/generax_56f3ed0:1.1.4'

    storeDir "${params.outdir}/species_trees/speciesrax"

    input:
    file input_bundles
    path reference_chronogram
    val species_names
    val outgroups

    output:
    path "inferred_species_tree.newick"  , emit: speciesrax_tree
    path "starting_species_tree.newick"
    path "mininj_species_tree.newick"
    path "species_tree_root_qc.tsv"
    path "species_tree_*.newick"
    path "*.txt"
    path "generax.log"
    path "mininj_generax.log"
    path "speciesrax_orthogroup.families"
    path "speciesrax_gene_tree_validation.tsv"
    path "speciesrax_selected_families.tsv"
    path "speciesrax_family_selection.tsv"
    path "speciesrax_selected_species_coverage.tsv"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (args.contains('--prune-species-tree') && args.contains('UndatedDL')) {
        throw new IllegalArgumentException(
            'GeneRax 2.1.3 cannot combine UndatedDL with --prune-species-tree: ' +
            'incomplete family coverage triggers an internal rate-vector assertion'
        )
    }
    """
    for archive in speciesrax_inputs_*.tar
    do
        tar -xf "\$archive"
        rm -f "\$archive"
    done

    prepare_speciesrax_gene_trees.py \\
        --manifest-glob 'speciesrax_inputs_*.tsv' \\
        --report speciesrax_gene_tree_validation.tsv

    echo "$species_names" \\
        | sed "s/\\[//g; s/\\]//g" \\
        | tr "," "\\n" \\
        | sed "s/^[[:space:]]*//; s/[[:space:]]*\$//" \\
        > speciesrax_expected_species.txt

    select_speciesrax_families.py \\
        --manifest-glob 'speciesrax_inputs_*.tsv' \\
        --validation-report speciesrax_gene_tree_validation.tsv \\
        --expected-species-file speciesrax_expected_species.txt \\
        --min-species-occupancy ${params.speciesrax_min_species_occupancy} \\
        --max-mean-copies ${params.speciesrax_max_mean_copies} \\
        --max-copies-per-species ${params.speciesrax_max_copies_per_species} \\
        --max-total-leaves-factor ${params.speciesrax_max_total_leaves_factor} \\
        --selected-manifest speciesrax_selected_families.tsv \\
        --report speciesrax_family_selection.tsv \\
        --species-coverage-report speciesrax_selected_species_coverage.tsv

    echo "[FAMILIES]" > speciesrax_orthogroup.families
    while IFS=\$'\\t' read -r og tree map_link
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

    # Pass 1: compute only the MPI-parallel MiniNJ topology. Its serialized
    # root is arbitrary, so do not optimize or use it for reconciliation.
    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree MiniNJ \\
        --families speciesrax_orthogroup.families \\
        --prefix MiniNJ \\
        --strategy SKIP \\
        --rec-model UndatedDL \\
        --si-strategy SKIP \\
        --do-not-reconcile

    cp MiniNJ/species_trees/inferred_species_tree.newick mininj_species_tree.newick
    cp MiniNJ/generax.log mininj_generax.log

    # Transfer a trusted root onto the exact MiniNJ topology. A reference
    # chronogram is preferred; explicit outgroups remain supported.
    if [ "${outgroups}" = "none" ]; then
        root_species_tree_from_reference.py \\
            --species-tree mininj_species_tree.newick \\
            --reference-chronogram ${reference_chronogram} \\
            --output rooted_mininj_species_tree.newick \\
            --qc-output species_tree_root_qc.tsv
    else
        root_species_tree_from_outgroups.py \\
            --species-tree mininj_species_tree.newick \\
            --outgroups '${outgroups}' \\
            --output rooted_mininj_species_tree.newick \\
            --qc-output species_tree_root_qc.tsv
    fi

    # Pass 2: preserve the rooted MiniNJ topology and root while estimating
    # final SpeciesRax branch lengths and quartet support under global D/L rates.
    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree rooted_mininj_species_tree.newick \\
        --families speciesrax_orthogroup.families \\
        --prefix SpeciesRax \\
        --strategy SKIP \\
        --si-estimate-bl \\
        $args

    mv SpeciesRax/* .
    rm -rf reconciliations results SpeciesRax MiniNJ
    mv species_trees/* .
    rm -r species_trees

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$(generax --version | head -n1 | sed 's/.*GeneRax //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
