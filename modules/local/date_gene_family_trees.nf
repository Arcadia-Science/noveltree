process DATE_GENE_FAMILY_TREES {
    tag "$meta.og"

    cpus 1
    // treePL only — PATHd8 is unsuitable for gene family trees (produces Inf
    // branches due to extreme rate heterogeneity between paralogs).
    // Retries with increasing time and memory; skip gracefully if all fail.
    time { 12.h * task.attempt }
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3
    memory {
        def n = (meta?.n_seq ?: 200) as long
        def estimated_gb = Math.max(4L, (long)(n * 4L / 1000L) + 1L)
        def capped_gb = (int) Math.min(estimated_gb, 24L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/zoogle:1.2.0'

    storeDir "${params.outdir}/gene_family_trees/time_calibrated"

    input:
    tuple val(meta), path(reconciled_tree), path(alignment), path(species_tree)
    val age_bracket

    output:
    tuple val(meta), path("${meta.og}_dated.newick")      , emit: dated_gft

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript /opt/zoogle/date_gene_family_tree.R \
        ${reconciled_tree} ${species_tree} ${alignment} \
        ${meta.og} ${age_bracket} \
        ${meta.og}_dated.newick ${meta.og}_calibrations.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -n1 | sed 's/R version //g' | cut -d' ' -f1 )
        ape: \$( Rscript -e "cat(as.character(packageVersion('ape')))" )
        phytools: \$( Rscript -e "cat(as.character(packageVersion('phytools')))" )
        phangorn: \$( Rscript -e "cat(as.character(packageVersion('phangorn')))" )
    END_VERSIONS
    """
}
