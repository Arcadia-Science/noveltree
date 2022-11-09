process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.3--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.3--hdfd78af_0' }"

    input:
    each mcl_inflation
    val ch_similarities

    output:
    path "*.fa" , emit: msas
    path "Orthogroups.GeneCount.tsv" , emit: genecounts
    path "Orthogroups.tsv" , emit: ogs_tsv
    path "Orthogroups.txt" , emit: ogs_txt
    path "Orthogroups_SingleCopyOrthologues.txt" , emit: scogs
    path "Orthogroups_UnassignedGenes.tsv" , emit: unassigned
    path "Orthogroups_SpeciesOverlaps.tsv" , emit: og_spp_overlap
    path "Statistics_Overall.tsv" , emit: stats
    path "Statistics_PerSpecies.tsv" , emit: stats_per_spp

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inflation_param = mcl_inflation ? "${mcl_inflation}" : '1.5'
    """
    orthofinder \\
        -b ../../../${params.outdir}/orthofinder/WorkingDirectory/ \\
        -n "Inflation_$inflation_param" \\
        -I $inflation_param \\
        -M msa -X -os -z \\
        -a $task.cpus \\
        $args
        
    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mv ../../../${params.outdir}/orthofinder/WorkingDirectory/OrthoFinder/Results_Inflation_$inflation_param/ ../../../${params.outdir}/orthofinder/
    
    # And clean up 
    rm -r ../../../${params.outdir}/orthofinder/WorkingDirectory/OrthoFinder/
    
    cp ../../../${params.outdir}/orthofinder/Results_Inflation_$inflation_param/*/*.* .
    """
}

