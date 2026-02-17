process FASTTREE {
    tag "$meta.og"
    label 'process_fasttree'

    container 'arcadiascience/fasttree_2.1.11:1.0.0'

    publishDir(
        path: "${params.outdir}/gene_family_trees/original",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), file(alignment)
    val model // not used

    output:
    tuple val(meta), path("*_ft.newick") , emit: phylogeny
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = alignment.baseName
    """
    # Make sure the number of threads are being specified properly
    export OMP_NUM_THREADS=${task.cpus}

    # Efficiently infer a gene family tree using FastTree!
    FastTreeDblMP \\
        $args \\
        $alignment > ${prefix}_ft.newick

    # prevent zero-length branches (sometimes inferred with fasttree)
    resolve_polytomies.R ${prefix}_ft.newick resolved.tree
    mv resolved.tree ${prefix}_ft.newick

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FastTree: \$(FastTreeDblMP 2>&1 | head -n1 | cut -d" " -f5)
    END_VERSIONS
    """
}
