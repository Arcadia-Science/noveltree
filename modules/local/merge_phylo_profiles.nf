process MERGE_PHYLO_PROFILES {
    tag "Merge Phylo Profiles"
    label "process_medium"

    container 'arcadiascience/phylo_profiles:1.0.0'

    input:
    path 'duplication_counts/*'
    path 'hgt_matrices/*'
    path 'loss_counts/*'
    path 'speciation_counts/*'
    path 'transfer_donor_counts/*'
    path 'transfer_recipient_counts/*'

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
    merge_phylo_profiles.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -n1 | sed "s/R version //g" | cut -f1 -d" " )
    END_VERSIONS
    """
}
