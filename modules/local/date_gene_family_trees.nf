process DATE_GENE_FAMILY_TREES {
    tag "$meta.og"
    label 'process_low_cpu'

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
