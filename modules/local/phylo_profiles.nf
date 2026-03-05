process PHYLO_PROFILES {
    tag "Phylo Profiles - ${meta.og}"
    label "process_medium"

    container 'arcadiascience/phylo_profiles:1.0.0'

    input:
    tuple val(meta), path(event_counts), path(species_event_counts), path(transfer_event_counts), path(species_coverage)
    path 'orthogroups'

    output:
    path "${meta.og}_duplication_count.tsv"        , emit: duplication_count
    path "${meta.og}_hgt_summed_counts.tsv"        , emit: hgt_summed_count
    path "${meta.og}_loss_count.tsv"               , emit: loss_count
    path "${meta.og}_speciation_count.tsv"         , emit: speciation_count
    path "${meta.og}_transfer_donor_count.tsv"     , emit: transfer_donor_count
    path "${meta.og}_transfer_recipient_count.tsv" , emit: transfer_recipient_count
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    phylo_profiles.R "${event_counts}" \\
                     "${species_event_counts}" \\
                     "${transfer_event_counts}" \\
                     "${species_coverage}" \\
                     "${meta.og}" \\
                     "orthogroups"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -n1 | sed "s/R version //g" | cut -f1 -d" " )
        data.table: \$( Rscript -e 'packageVersion("data.table")' | cut -f2 -d" " | sed "s/'//g" | sed "s/'//g" )
        plyr: \$( Rscript -e 'packageVersion("plyr")' | cut -f2 -d" " | sed "s/'//g" | sed "s/'//g" )
        purrr: \$( Rscript -e 'packageVersion("purrr")' | cut -f2 -d" " | sed "s/'//g" | sed "s/'//g" )
    END_VERSIONS
    """
}
