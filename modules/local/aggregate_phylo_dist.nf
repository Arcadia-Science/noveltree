process AGGREGATE_PHYLO_DIST {
    tag "Aggregate PHYLO_DIST results"
    label "process_low"

    publishDir(
        path: "${params.outdir}/phylo_dist_aggregated",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    )

    input:
    path summary_tables  // All OG*_final_summary_table.tsv files

    output:
    path "all_protein_comparisons.tsv.gz", emit: aggregated_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Count input files for logging
    echo "Aggregating \$(ls -1 *.tsv | wc -l) gene family tables..."

    # Concatenate all tables with single header, compress
    (
        # Take header from first file
        head -1 \$(ls *.tsv | head -1)
        # Append all data lines (skip headers with tail -n +2)
        tail -n +2 -q *.tsv
    ) | gzip -c > all_protein_comparisons.tsv.gz

    echo "Aggregation complete. Output: all_protein_comparisons.tsv.gz"

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(gzip --version 2>&1 | head -1 | sed 's/gzip //' | cut -d' ' -f1)
    END_VERSIONS
    """
}