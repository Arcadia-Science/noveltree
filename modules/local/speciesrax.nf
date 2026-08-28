process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_generax'
    stageInMode 'copy'

    // 1.1.4 adds Biopython for deterministic in-process root transfer.
    container 'arcadiascience/generax_56f3ed0:1.1.4'

    storeDir "${params.outdir}/species_trees/speciesrax"

    input:
    file input_bundles
    // Nextflow represents a missing optional file as an empty list. The
    // outgroup branch passes [] and never dereferences this input.
    path reference_chronogram
    val species_names
    val outgroups
    path speciesrax_selection
    path speciesrax_coverage

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

    # Family eligibility is already determined from the native OrthoFinder
    # gene-count table. Consolidate the bounded bundle manifests without
    # repeating selection or repairing task-local mapping files.
    printf 'orthogroup\tgene_tree\tmapping\n' > speciesrax_selected_families.tsv
    for manifest in \$(find . -maxdepth 1 -name 'speciesrax_inputs_*.tsv' -print | LC_ALL=C sort)
    do
        tail -n +2 "\$manifest"
    done | LC_ALL=C sort -k1,1 >> speciesrax_selected_families.tsv

    duplicate_ogs=\$(tail -n +2 speciesrax_selected_families.tsv | cut -f1 | uniq -d | head)
    if [ -n "\$duplicate_ogs" ]; then
        echo "Duplicate SpeciesRax family manifests: \$duplicate_ogs" >&2
        exit 1
    fi
    if [ "\$(wc -l < speciesrax_selected_families.tsv)" -le 1 ]; then
        echo "No upstream-selected SpeciesRax families were staged" >&2
        exit 1
    fi
    # These audit tables are already staged under their final output names.
    # Validate them in place instead of copying a file onto itself, which GNU
    # cp treats as an error under `set -e`.
    for upstream_audit in "${speciesrax_selection}" "${speciesrax_coverage}"
    do
        if [ ! -s "\$upstream_audit" ]; then
            echo "Missing or empty upstream SpeciesRax audit: \$upstream_audit" >&2
            exit 1
        fi
    done

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
    rm -f speciesrax_inputs_*.tsv

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
    validate_root_split.py \
        --expected rooted_mininj_species_tree.newick \
        --observed inferred_species_tree.newick \
        --qc-output species_tree_root_qc.tsv
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
