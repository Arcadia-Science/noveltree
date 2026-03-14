process TIME_CALIBRATE_SPECIES_TREE {
    tag "Time-calibrate species tree"
    label 'process_single'

    container 'arcadiascience/zoogle:1.0.0'

    storeDir "${params.outdir}/species_trees/time_calibrated"

    input:
    path species_tree       // Species tree from SpeciesRax (Newick format)
    path reference_tree     // User-provided time-calibrated reference tree (Newick format)
    val calibration_method  // Method for time calibration: "treePL" or "PATHd8"
    val age_bracket         // Fractional uncertainty for calibration ages (e.g. 0.20)

    output:
    path "time_calibrated_species_tree.newick", emit: calibrated_tree
    path "calibration_log.txt"                , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Run the R script to time-calibrate the species tree
    # The R script normalizes reference tree tips to hyphens internally,
    # so the output already uses the pipeline's hyphen convention.
    time_calibrate_species_tree.R \\
        ${species_tree} \\
        ${reference_tree} \\
        time_calibrated_species_tree.newick \\
        ${calibration_method} \\
        ${age_bracket} \\
        2>&1 | tee calibration_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //g' | cut -d' ' -f1)
        ape: \$(Rscript -e "cat(as.character(packageVersion('ape')))")
        phangorn: \$(Rscript -e "cat(as.character(packageVersion('phangorn')))")
    END_VERSIONS
    """
}
