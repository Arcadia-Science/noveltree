process FASTTREE {
    tag "$alignment"
    label 'process_iqtree'

    container "${ workflow.containerEngine == 'docker' ? 'staphb/fasttree:2.1.11':
        '' }"

    publishDir(
        path: "${params.outdir}/fasttree_gene_trees",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
    input:
    path(alignment)
    val model // not used

    output:
    path("*.treefile")  , emit: phylogeny
    path(alignment)     , emit: msa
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    """
    og=\$(echo $alignment | cut -f1 -d "_")
    
    # Efficiently infer a gene family tree using FastTree!
    /MAGUS/magus_tools/fasttree/FastTreeMP \\
        $args \\
        $alignment > \${og}_ft.treefile
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FastTree: \$(/MAGUS/magus_tools/fasttree/FastTreeMP | head -n1 | cut -d" " -f5)
    END_VERSIONS
    """
}
