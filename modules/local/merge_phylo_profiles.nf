// Sum HGT matrices element-wise across per-OG outputs.
// TSV concatenation is handled upstream via collectFile() to avoid
// staging hundreds of small files from S3 into a single process.
process MERGE_PHYLO_PROFILES {
    tag "Sum HGT matrices"
    label "process_single"

    container 'arcadiascience/phylo_profiles:1.0.0'

    input:
    path 'hgt_matrices/*'

    output:
    path "hgt_summed_counts_recip_donor.tsv" , emit: hgt_summed_count
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sum_hgt_matrices.py hgt_matrices hgt_summed_counts_recip_donor.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | cut -d' ' -f2 )
    END_VERSIONS
    """
}
