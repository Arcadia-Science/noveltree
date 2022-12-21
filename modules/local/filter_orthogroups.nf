process FILTER_ORTHOGROUPS {
    tag "Summarize OG taxon distribution"
    label 'process_single'

    // Dockerhubs r-base container doesn't have ps (procps)
    // installed, which is required by nextflow to monitor
    // processes. So, we use the cogeqc module here for simplicity.
    container "${ workflow.containerEngine == 'docker' ?
        'austinhpatton123/cogeqc-1.2.0_r-4.2.2': '' }"

    publishDir(
        path: "${params.outdir}/orthofinder/full_analysis",
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
    path "all_ogs_counts.csv"           , emit: all_ogs
    path "spptree_core_ogs_counts.csv"  , emit: spptree_core_ogs
    path "genetree_core_ogs_counts.csv" , emit: genetree_core_ogs
    path "SpeciesTreeFastas/*.fa"       , emit: spptree_fas
    path "GeneTreeFastas/*.fa"          , emit: genetree_fas

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def num_spp = min_num_spp ? "${min_num_spp}" : '4'
    def num_grp = min_num_groups ? "${min_num_groups}" : '2'
    def copy_num1 = max_copy_num_filt1 ? "${max_copy_num_filt1}" : '5'
    def copy_num2 = max_copy_num_filt2 ? "${max_copy_num_filt2}" : '10'

    """
    # The below summarizes the distribution of orthogroups across
    # species and across taxonomic groups, filtering them for
    # species/gene tree inference based on their distribution
    # and copy numbers

    ogSppCounts=${orthofinder_outdir}/Orthogroups/Orthogroups.GeneCount.tsv

    # Run the scripts to generate the orthogroup species/taxa gene count summaries and filtered sets
    Rscript $projectDir/bin/og-tax-summary.R \$ogSppCounts $samplesheet $num_spp $num_grp $copy_num1 $copy_num2

    # Add a column of filepaths to these
    msaDir=\$( cd ${orthofinder_outdir}/Orthogroup_Sequences/; pwd )

    # Create the column
    tail -n+2 spptree_core_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msaDir}/&.fa|g" | sed '1 i\\file' > spptree_core_og_fpaths.txt
    tail -n+2 genetree_core_ogs_counts.csv | cut -f1 -d"," | sed "s|.*|\${msaDir}/&.fa|g" | sed '1 i\\file' > genetree_core_og_fpaths.txt

    mkdir SpeciesTreeFastas
    mkdir GeneTreeFastas

    while IFS= read -r line
    do
        # Copy the file to the destination directory
        cp "\$line" SpeciesTreeFastas/
    done < spptree_core_og_fpaths.txt
    while IFS= read -r line
    do
        # Copy the file to the destination directory
        cp "\$line" GeneTreeFastas/
    done < genetree_core_og_fpaths.txt
    """
}
