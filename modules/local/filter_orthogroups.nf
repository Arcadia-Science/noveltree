process FILTER_ORTHOGROUPS {
    tag "Summarize OG taxon distribution"
    label 'process_single'

    // Dockerhubs r-base container doesn't have ps (procps)
    // installed, which is required by nextflow to monitor
    // processes. So, we use the cogeqc module here for simplicity.
    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/rbase:4.2.2': '' }"

    publishDir(
        path: "${params.outdir}/filtered_orthogroups",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path samplesheet        // Path to samplesheet produced by input check containing sample metadata
    path orthofinder_outdir // Directory containing all the inflation params
    val min_num_spp         // Minimum number of species to infer MSAs/trees for
    val min_num_groups      // Minimum number of clades/taxonomic groups
    val max_copy_num_filt1  // Max copy number for genes intended for species tree inference
    val max_copy_num_filt2  // Max copy number for all other gene family tree inference

    output:
    path "*"                            , emit: filtered_ogs
    path "all_ogs_counts.csv"           , emit: all_ogs
    path "spptree_core_ogs_counts.csv"  , emit: spptree_core_ogs
    path "genetree_core_ogs_counts.csv" , emit: genetree_core_ogs
    path "species_tree_og_msas/*.fa"    , emit: spptree_fas
    path "gene_tree_og_msas/*.fa"       , emit: genetree_fas

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # The below summarizes the distribution of orthogroups across
    # species and across taxonomic groups, filtering them for
    # species/gene tree inference based on their distribution
    # and copy numbers

    og_spp_counts=${orthofinder_outdir}/Orthogroups/Orthogroups.GeneCount.tsv

    # Run the scripts to generate the orthogroup species/taxa gene count summaries and filtered sets
    og_tax_summary.R \$og_spp_counts $samplesheet $min_num_spp $min_num_groups $max_copy_num_filt1 $max_copy_num_filt2

    # Add a column of filepaths to these
    msa_dir=\$( cd ${orthofinder_outdir}/Orthogroup_Sequences/; pwd )

    # Create the column
    tail -n+2 spptree_core_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msa_dir}/&.fa|g" > spptree_core_og_fpaths.txt
    tail -n+2 genetree_core_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msa_dir}/&.fa|g" > genetree_core_og_fpaths.txt

    mkdir species_tree_og_msas
    mkdir gene_tree_og_msas

    while IFS= read -r trees
    do
        # Copy the file to the destination directory
        cp "\$trees" species_tree_og_msas/
    done < spptree_core_og_fpaths.txt
    while IFS= read -r trees
    do
        # Copy the file to the destination directory
        cp "\$trees" gene_tree_og_msas/
    done < genetree_core_og_fpaths.txt

    # Remove these intermediate files
    rm *fpaths.txt
    """
}
