process PHYLO_PROFILES {
    tag "Phylo Profiles"
    label "process_high_memory"

    container 'arcadiascience/phylo_profiles:1.0.0'

    input:
    path 'event_count_files'
    path 'species_event_counts_files'
    path 'transfer_event_counts_files'
    path 'species_coverage_files'
    val 'ogs'
    path 'orthogroups'

    output:
    path "duplication_count_per_species_per_gene_family.tsv"        , emit: duplication_count
    path "hgt_summed_counts_recip_donor.tsv"                        , emit: hgt_summed_count
    path "loss_count_per_per_species_gene_family.tsv"               , emit: loss_count
    path "speciation_count_per_species_per_gene_family.tsv"         , emit: speciation_count
    path "transfer_donor_count_per_species_per_gene_family.tsv"     , emit: transfer_donor_count
    path "transfer_recipient_count_per_species_per_gene_family.tsv" , emit: transfer_recipient_count
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    phylo_profiles.R "${event_count_files}" \\
                     "${species_event_counts_files}" \\
                     "${transfer_event_counts_files}" \\
                     "${species_coverage_files}" \\
                     "${ogs.join(' ')}" \\
                     "${orthogroups}"

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
