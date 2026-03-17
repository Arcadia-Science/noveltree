process DATE_GENE_FAMILY_TREES {
    tag "$meta.og"

    cpus 1
    // Attempt 1: treePL, 12h, base memory
    // Attempt 2: PATHd8, 12h, same memory
    // Attempt 3: PATHd8, 24h, 2x memory
    // Attempt 4: PATHd8, 24h, 3x memory
    time { task.attempt <= 2 ? 12.h : 24.h }
    errorStrategy { task.attempt <= 4 ? 'retry' : 'ignore' }
    maxRetries 4
    memory {
        def n = (meta?.n_seq ?: 200) as long
        def estimated_gb = Math.max(2L, (long)(n * 4L / 1000L) + 1L)
        def capped_gb = (int) Math.min(estimated_gb, 12L)
        def mem_multiplier = task.attempt <= 2 ? 1 : (task.attempt - 1)
        def requested = capped_gb.GB * mem_multiplier
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/zoogle:1.0.0'

    storeDir "${params.outdir}/gene_family_trees/time_calibrated"

    input:
    tuple val(meta), path(reconciled_tree), path(alignment), path(species_tree)
    val max_treepl_tips

    output:
    tuple val(meta), path("${meta.og}_dated.newick")      , emit: dated_gft

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // On retry, force PATHd8 by setting max_treepl_tips to 0
    def effective_max_tips = task.attempt > 1 ? 0 : max_treepl_tips
    """
    Rscript /opt/zoogle/date_gene_family_tree.R \
        ${reconciled_tree} ${species_tree} ${alignment} \
        ${meta.og} ${effective_max_tips} \
        ${meta.og}_dated.newick ${meta.og}_calibrations.csv
    """
}
