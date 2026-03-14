process ZOOGLE {
    tag "${meta.og}"
    label "process_high"

    container 'arcadiascience/zoogle:1.0.0'

    storeDir "${params.outdir}/zoogle"

    input:
    tuple val(meta), path(gene_tree), path(phys_props_file)
    val ref_species

    output:
    path "phylo-corrected-data/${meta.og}_phylo_corr_dat.tsv"           , emit: phylo_corrected_data
    path "protein-dist-mats/${meta.og}_protein_dists.tsv"               , emit: protein_dist_mat
    path "protein-phylo-dist-mats/${meta.og}_phylo_dists.tsv"           , emit: prot_phylo_dists
    path "protein-dists-to-reference/${meta.og}_protein_dists.tsv"      , emit: prot_dists_to_ref
    path "species-dists-to-reference/${meta.og}_species_dists.tsv"      , emit: spp_dists_to_ref
    path "protein-pvals/${meta.og}_protein_reference_dist_pvals.tsv"    , emit: protein_pvals
    path "species-pvals/${meta.og}_species_reference_dist_pvals.tsv"    , emit: species_pvals
    path "pairwise-protein-dist-perm-test/${meta.og}_protein_protein_dist_permutation_test.tsv" , emit: per_protein_dist_res
    path "final_protein_pair_summary_tables/${meta.og}_final_summary_table.tsv" , emit: final_summary_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript

    # Enable error tracing
    options(error = function() {
        traceback(2)
        quit(status = 1)
    })

    # Save the work directory path where input files are staged
    work_dir <- getwd()
    cat("Work directory:", work_dir, "\\n")

    # Change to /opt/zoogle to source R scripts (they need to source C++ files with relative paths)
    setwd("/opt/zoogle")
    cat("Sourcing R scripts from /opt/zoogle...\\n")
    source("phylo_multivariate_distance_functions.R")
    source("protein_dist_permutation_tests.R")
    source("simple_protein_dist_signif_tests.R")
    source("protein_distance_calculation_functions.R")
    cat("Scripts sourced successfully\\n")

    # Change back to work directory where input files are
    setwd(work_dir)
    cat("Changed back to work directory\\n")

    # Prepare gene family info as named vector
    gene_family <- c(
        family = "${meta.og}",
        gft = "${gene_tree}"
    )

    # Call the main function
    # Note: gene_family["family"] is the OG name, script expects aa_stat_basedir/OG_summary_statistics.csv
    # We need to point to current dir since Nextflow stages the file here
    genefam_aa_conservation(
        gene_family = gene_family,
        ref_spp = "${ref_species}",
        aa_stat_basedir = "",
        keep_stats = c("molecular_weight", "aromaticity", "instability", "flexibility",
                       "gravy_bm", "isoelectric_point", "charge_at_pH_7", "helix_fract",
                       "sheet_fract", "molar_ext_coef_cysteines"),
        out_dir = "."
    )

    # Create versions file
    writeLines(
        c(
            '"${task.process}":',
            paste0('    R: "', R.version.string, '"'),
            paste0('    Rcpp: "', packageVersion("Rcpp"), '"'),
            paste0('    RcppArmadillo: "', packageVersion("RcppArmadillo"), '"'),
            paste0('    ape: "', packageVersion("ape"), '"'),
            paste0('    phytools: "', packageVersion("phytools"), '"'),
            paste0('    geiger: "', packageVersion("geiger"), '"')
        ),
        "versions.yml"
    )
    """
}
