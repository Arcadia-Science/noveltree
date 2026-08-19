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

    output:
    path("*/Results_Inflation*"),           emit: inflation_dir
    path("species_tree_og_fas/*.fa"),      emit: spptree_fas, optional: true
    path("gene_tree_og_fas/*.fa"),         emit: genetree_fas, optional: true
    path("all_ogs_counts.csv"),             emit: all_ogs, optional: true
    path("spptree_core_ogs_counts.csv"),    emit: spptree_core_ogs, optional: true
    path("genetree_core_ogs_counts.csv"),   emit: genetree_core_ogs, optional: true
    path("og_fasta_metadata.tsv"),          emit: og_metadata, optional: true
    path("chimera_report.tsv"),             emit: chimera_report, optional: true
    path("orthofinder_checkpoint_receipt.json"), emit: checkpoint_receipt, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def checkpoint_root = "${params.outdir}/orthofinder/mcl_checkpoints"
    """
    checkpoint_fingerprint=''
    checkpoint_restored=false
    shopt -s nullglob

    # Hash the biological/clustering inputs before deleting the bundles. This
    # key intentionally ignores NovelTree postprocessing code so a retry after
    # a postprocessing fix can reuse the completed OrthoFinder result.
    if [ "$output_directory" == "complete_dataset" ]; then
        checkpoint_inputs=(SpeciesIDs.txt SequenceIDs.txt Species*.fa)
        checkpoint_bundles=(BlastBundle_*.tar)
        if [ "\${#checkpoint_inputs[@]}" -lt 3 ] || [ "\${#checkpoint_bundles[@]}" -lt 1 ]; then
            echo "Insufficient inputs to fingerprint the OrthoFinder clustering" >&2
            exit 1
        fi
        checkpoint_fingerprint=\$(orthofinder_checkpoint.py fingerprint \
            --orthofinder-version '2.5.4' \
            --inflation '${mcl_inflation}' \
            --orthofinder-options '-b -I -M msa -X -os -z' \
            --extra-args '${args}' \
            "\${checkpoint_inputs[@]}" \
            --metadata-only-inputs "\${checkpoint_bundles[@]}")
    fi

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

    # Restore raw tables and orthogroup FASTAs when an earlier attempt reached
    # the end of OrthoFinder but failed during custom postprocessing.
    if [ "$output_directory" == "complete_dataset" ]; then
        if orthofinder_checkpoint.py restore \
            --root '${checkpoint_root}' \
            --fingerprint "\$checkpoint_fingerprint" \
            --results-dir 'OrthoFinder/Results_Inflation_${mcl_inflation}' \
            --receipt orthofinder_checkpoint_receipt.json
        then
            checkpoint_restored=true
        else
            restore_status=\$?
            if [ "\$restore_status" -ne 10 ]; then
                exit "\$restore_status"
            fi
        fi
    fi

    if [ "\$checkpoint_restored" != true ]; then
        orthofinder \\
            -b ./ \\
            -n "Inflation_${mcl_inflation}" \\
            -I $mcl_inflation \\
            -M msa -X -os -z \\
            -t ${task.cpus} \\
            -a ${task.cpus} \\
            $args

        # Persist the unmodified clustering immediately after OrthoFinder
        # succeeds, before any custom filtering can fail or mutate it.
        if [ "$output_directory" == "complete_dataset" ]; then
            orthofinder_checkpoint.py save \
                --root '${checkpoint_root}' \
                --fingerprint "\$checkpoint_fingerprint" \
                --results-dir 'OrthoFinder/Results_Inflation_${mcl_inflation}' \
                --receipt orthofinder_checkpoint_receipt.json
        fi
    fi

    # Check if we're running an mcl test or not:
    # if so, delete the sequence files and other non-essential directories that
    # we will not be using and take up significant, unnecessary space.
    if [ "$output_directory" == "mcl_test_dataset" ]; then
        rm -r OrthoFinder/*/Single_Copy_Orthologue_Sequences/
        rm -r OrthoFinder/*/Orthogroup_Sequences/
        rm -r OrthoFinder/*/WorkingDirectory/
        rm -r OrthoFinder/*/Orthologues/
    else
        # Preliminary streaming filter. Chimera detection only needs to score
        # proteins in OGs that can proceed downstream.
        og_tax_summary.py \\
            OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.GeneCount.tsv \\
            ${samplesheet} \\
            ${min_num_seqs} ${min_num_spp} ${min_prop_spp_for_spptree} ${max_copy_num}

        # Flag and remove cross-OG chimeric proteins. All protein assignments
        # remain available as possible hit OGs, but queries and FASTA rewrites
        # are restricted to preliminary retained families.
        flag_cross_og_chimeras.py \\
            --orthogroups OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.tsv \\
            --blast_dir ./ \\
            --og_seqs_dir OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroup_Sequences/ \\
            --report chimera_report.tsv \\
            --retained-og-files spptree_core_ogs_counts.csv genetree_core_ogs_counts.csv \\
            --gene-counts OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.GeneCount.tsv \\
            --updated-gene-counts Orthogroups.GeneCount.post_chimera.tsv \\
            --updated-orthogroups Orthogroups.post_chimera.tsv

        # Recompute the final retained sets after chimera removal. An excluded
        # OG cannot become eligible when sequences are removed, so the
        # preliminary restriction above is lossless.
        og_tax_summary.py \\
            Orthogroups.GeneCount.post_chimera.tsv \\
            ${samplesheet} \\
            ${min_num_seqs} ${min_num_spp} ${min_prop_spp_for_spptree} ${max_copy_num}

        # Ensure the persisted OrthoFinder directory and downstream
        # PHYLO_PROFILES use membership/count tables matching the filtered
        # FASTAs rather than the pre-chimera OrthoFinder tables.
        mv Orthogroups.post_chimera.tsv \\
            OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.tsv
        mv Orthogroups.GeneCount.post_chimera.tsv \\
            OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.GeneCount.tsv

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
