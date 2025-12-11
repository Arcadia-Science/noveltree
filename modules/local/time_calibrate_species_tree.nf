process TIME_CALIBRATE_SPECIES_TREE {
    tag "Time-calibrate species tree"
    label 'process_single'

    container 'arcadiascience/phylo_dist:1.0.0'

    publishDir(
        path: "${params.outdir}/time_calibrated_species_tree",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path species_tree       // Species tree from SpeciesRax (Newick format)
    path reference_tree     // User-provided time-calibrated reference tree (Newick format)
    val calibration_method  // Method for time calibration: "treePL" or "PATHd8"

    output:
    path "time_calibrated_species_tree.newick", emit: calibrated_tree
    path "calibration_log.txt"                , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Run the R script to time-calibrate the species tree
    time_calibrate_species_tree.R \\
        ${species_tree} \\
        ${reference_tree} \\
        time_calibrated_species_tree.newick \\
        ${calibration_method} \\
        2>&1 | tee calibration_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //g' | cut -d' ' -f1)
        ape: \$(Rscript -e "cat(as.character(packageVersion('ape')))")
        geiger: \$(Rscript -e "cat(as.character(packageVersion('geiger')))")
    END_VERSIONS
    """
}
