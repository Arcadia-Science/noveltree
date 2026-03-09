process DATE_GENE_FAMILY_TREES {
    tag "$meta.og"

    cpus 1
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 200) as long
        def estimated_gb = Math.max(2L, (long)(n * 4L / 1000L) + 1L)
        def capped_gb = (int) Math.min(estimated_gb, 12L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/zoogle:1.0.0'

    input:
    tuple val(meta), path(reconciled_tree), path(alignment), path(species_tree)
    val max_treepl_tips

    output:
    tuple val(meta), path("${meta.og}_dated.newick")      , emit: dated_gft
    tuple val(meta), path("${meta.og}_calibrations.csv")   , emit: calibrations
    path "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript /opt/zoogle/date_gene_family_tree.R \
        ${reconciled_tree} ${species_tree} ${alignment} \
        ${meta.og} ${max_treepl_tips} \
        ${meta.og}_dated.newick ${meta.og}_calibrations.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //g' | cut -d' ' -f1)
        ape: \$(Rscript -e "cat(as.character(packageVersion('ape')))")
        phangorn: \$(Rscript -e "cat(as.character(packageVersion('phangorn')))")
        phytools: \$(Rscript -e "cat(as.character(packageVersion('phytools')))")
    END_VERSIONS
    """
}
