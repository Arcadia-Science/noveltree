// Pivot the concatenated long-format HGT counts into a species × species matrix.
// All other phylo profile merges are handled by collectFile() in the workflow.
process MERGE_PHYLO_PROFILES {
    tag "Pivot HGT matrix"
    label "process_single"

    container 'arcadiascience/phylo_profiles:1.0.0'

    storeDir "${params.outdir}/gene_family_evolution"

    input:
    path hgt_long_file

    output:
    path "hgt_summed_counts_recip_donor.tsv" , emit: hgt_summed_count

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pivot_hgt_matrix.py ${hgt_long_file} hgt_summed_counts_recip_donor.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | cut -d' ' -f2 )
    END_VERSIONS
    """
}
