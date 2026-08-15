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
    path("chimera_report.tsv"),             emit: chimera_report, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def analysis_cpus = Math.max(1, (task.cpus as int).intdiv(4))

    """
    # Expand one archive per query species, removing each archive immediately
    # to avoid retaining both the archive and extracted files on local scratch.
    for archive in BlastBundle_*.tar
    do
        tar -xf \$archive
        rm -f \$archive
    done

    for f in \$(ls TestBlast*)
    do
        mv \$f \$(echo \$f | sed "s/TestBlast/Blast/g")
    done

    orthofinder \\
        -b ./ \\
        -n "Inflation_${mcl_inflation}" \\
        -I $mcl_inflation \\
        -M msa -X -os -z \\
        -t ${task.cpus} \\
        -a ${analysis_cpus} \\
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
        dir=\$(pwd)
        cd \$(ls -d OrthoFinder/*/WorkingDirectory)
        tar -czvf Sequences_ids.tar.gz Sequences_ids
        rm -r Sequences_ids
        cd \$dir

        # Flag and remove cross-OG chimeric proteins
        flag_cross_og_chimeras.py \\
            --orthogroups OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.tsv \\
            --blast_dir ./ \\
            --og_seqs_dir OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroup_Sequences/ \\
            --report chimera_report.tsv

        # Filter orthogroups into species tree and gene tree sets
        og_tax_summary.py \\
            OrthoFinder/Results_Inflation_${mcl_inflation}/Orthogroups/Orthogroups.GeneCount.tsv \\
            ${samplesheet} \\
            ${min_num_seqs} ${min_num_spp} ${min_prop_spp_for_spptree} ${max_copy_num}

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
