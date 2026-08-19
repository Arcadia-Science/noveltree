process PHYLO_PROFILES {
    tag "Phylo Profiles - ${meta.og}"
    label "process_single"

    container 'arcadiascience/phylo_profiles:1.0.0'

    input:
    tuple val(meta), path(species_event_counts), path(species_coverage)

    output:
    path "${meta.og}_duplication_count.tsv" , emit: duplication_count
    path "${meta.og}_loss_count.tsv"        , emit: loss_count
    path "${meta.og}_speciation_count.tsv" , emit: speciation_count

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    phylo_profiles.R "${species_event_counts}" "${species_coverage}" "${meta.og}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -n1 | sed "s/R version //g" | cut -f1 -d" " )
    END_VERSIONS
    """
}
