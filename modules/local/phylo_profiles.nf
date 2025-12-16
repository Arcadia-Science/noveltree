process PHYLO_PROFILES {
    tag "Phylo Profiles - Batch ${batch_id}"
    label "process_medium"

    container 'arcadiascience/phylo_profiles:1.0.0'

    input:
    tuple val(batch_id), path(event_count_files), path(species_event_counts_files), path(transfer_event_counts_files), path(species_coverage_files), val(ogs)
    path 'orthogroups'

    output:
    path "batch_${batch_id}_duplication_count.tsv"        , emit: duplication_count
    path "batch_${batch_id}_hgt_summed_counts.tsv"        , emit: hgt_summed_count
    path "batch_${batch_id}_loss_count.tsv"               , emit: loss_count
    path "batch_${batch_id}_speciation_count.tsv"         , emit: speciation_count
    path "batch_${batch_id}_transfer_donor_count.tsv"     , emit: transfer_donor_count
    path "batch_${batch_id}_transfer_recipient_count.tsv" , emit: transfer_recipient_count
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def event_files = event_count_files.collect { it.name }.join(' ')
    def species_event_files = species_event_counts_files.collect { it.name }.join(' ')
    def transfer_files = transfer_event_counts_files.collect { it.name }.join(' ')
    def coverage_files = species_coverage_files.collect { it.name }.join(' ')
    def og_list = ogs.join(' ')
    """
    phylo_profiles.R "${event_files}" \\
                     "${species_event_files}" \\
                     "${transfer_files}" \\
                     "${coverage_files}" \\
                     "${og_list}" \\
                     "orthogroups" \\
                     "${batch_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -n1 | sed "s/R version //g" | cut -f1 -d" " )
        data.table: \$( Rscript -e 'packageVersion("data.table")' | cut -f2 -d" " | sed "s/‘//g" | sed "s/’//g" )
        parallel: \$( Rscript -e 'packageVersion("parallel")' | cut -f2 -d" " | sed "s/‘//g" | sed "s/’//g" )
        plyr: \$( Rscript -e 'packageVersion("plyr")' | cut -f2 -d" " | sed "s/‘//g" | sed "s/’//g" )
        purrr: \$( Rscript -e 'packageVersion("purrr")' | cut -f2 -d" " | sed "s/‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}
