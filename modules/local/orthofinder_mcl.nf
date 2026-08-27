process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_orthofinder'

    storeDir "${params.outdir}/orthofinder/mcl"

    container 'arcadiascience/orthofinder_2.5.4:1.0.0'

    input:
    each mcl_inflation
    file(blast_bundles)
    file(fasta)
    file(db)
    file(sppIDs)
    file(seqIDs)
    val output_directory
    file samplesheet
    val min_num_seqs
    val min_num_spp
    val min_prop_spp_for_spptree
    val max_copy_num
    path canonical_maps, arity: '0..*'

    output:
    path("*/Results_Inflation*"),           emit: inflation_dir
    path("species_tree_og_fas/*.fa"),      emit: spptree_fas, optional: true
    path("gene_tree_og_fas/*.fa"),         emit: genetree_fas, optional: true
    path("all_ogs_counts.csv"),             emit: all_ogs, optional: true
    path("spptree_core_ogs_counts.csv"),    emit: spptree_core_ogs, optional: true
    path("genetree_core_ogs_counts.csv"),   emit: genetree_core_ogs, optional: true
    path("og_fasta_metadata.tsv"),          emit: og_metadata, optional: true
    path("family_maps/*.map.link"),         emit: family_maps, optional: true
    path("speciesrax_family_selection.tsv"), emit: speciesrax_selection, optional: true
    path("speciesrax_selected_species_coverage.tsv"), emit: speciesrax_coverage, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    shopt -s nullglob

    # Expand one archive per query species, removing each archive immediately
    # to avoid retaining both the archive and extracted files on local scratch.
    for archive in BlastBundle_*.tar
    do
        tar -xf \$archive
        rm -f \$archive
    done

    for f in TestBlast*
    do
        mv \$f \$(echo \$f | sed "s/TestBlast/Blast/g")
    done

    orthofinder \\
        -b ./ \\
        -n "Inflation_${mcl_inflation}" \\
        -I $mcl_inflation \\
        -M msa -X -os -z \\
        -t ${task.cpus} \\
        -a ${task.cpus} \\
        $args

    # Check if we're running an mcl test or not:
    # if so, delete the sequence files and other non-essential directories that
    # we will not be using and take up significant, unnecessary space.
    if [ "$output_directory" == "mcl_test_dataset" ]; then
        rm -r OrthoFinder/*/Single_Copy_Orthologue_Sequences/
        rm -r OrthoFinder/*/Orthogroup_Sequences/
        rm -r OrthoFinder/*/WorkingDirectory/
        rm -r OrthoFinder/*/Orthologues/
    else
        # Filter the native OrthoFinder gene-count table once. Orthogroup
        # membership is not mutated after clustering.
        og_tax_summary.py \\
            OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.GeneCount.tsv \\
            ${samplesheet} \\
            ${min_num_seqs} ${min_num_spp} ${min_prop_spp_for_spptree} ${max_copy_num} \\
            ${params.speciesrax_min_species_occupancy} \\
            ${params.speciesrax_max_mean_copies} \\
            ${params.speciesrax_max_copies_per_species} \\
            ${params.speciesrax_max_total_leaves_factor}

        # Move filtered FASTAs into separate directories
        msa_dir=OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroup_Sequences

        mkdir -p species_tree_og_fas gene_tree_og_fas

        tail -n+2 spptree_core_ogs_counts.csv | cut -f1 -d"," | while read og; do
            if [ -f "\${msa_dir}/\${og}.fa" ]; then
                mv "\${msa_dir}/\${og}.fa" species_tree_og_fas/
            fi
        done

        tail -n+2 genetree_core_ogs_counts.csv | cut -f1 -d"," | while read og; do
            if [ -f "\${msa_dir}/\${og}.fa" ]; then
                mv "\${msa_dir}/\${og}.fa" gene_tree_og_fas/
            fi
        done

        # Join retained proteins to the canonical input mapping once. Every
        # downstream aligner and trimmer propagates or subsets these files.
        build_family_maps.py \
            --canonical-map-glob '*.protein_map.tsv' \
            --species-tree-dir species_tree_og_fas \
            --gene-tree-dir gene_tree_og_fas \
            --output-dir family_maps

        # Summarize the retained FASTAs once on local scratch. Downstream
        # resource routing consumes this small manifest instead of asking the
        # Nextflow controller to download and parse every FASTA from S3.
        summarize_og_fastas.py \
            --species-tree-dir species_tree_og_fas \
            --gene-tree-dir gene_tree_og_fas \
            --output og_fasta_metadata.tsv

        # Remove directories no longer needed (orthology derived from GeneRax reconciliations in PARSE_PHYLOHOGS)
        rm -rf OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroup_Sequences/
        rm -rf OrthoFinder/Results_Inflation_${mcl_inflation}/Single_Copy_Orthologue_Sequences/
        rm -rf OrthoFinder/Results_Inflation_${mcl_inflation}/Orthologues/
        rm -rf OrthoFinder/Results_Inflation_${mcl_inflation}/Comparative_Genomics_Statistics/
        rm -rf OrthoFinder/Results_Inflation_${mcl_inflation}/Gene_Trees/
        rm -rf OrthoFinder/Results_Inflation_${mcl_inflation}/WorkingDirectory/
        rm -f  OrthoFinder/Results_Inflation_${mcl_inflation}/Citation.txt
    fi

    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mkdir ${output_directory}
    mv OrthoFinder/Results_Inflation_${mcl_inflation}/ ${output_directory}/
    rm -r OrthoFinder/
    """
}
