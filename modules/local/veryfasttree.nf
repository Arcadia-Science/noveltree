process VERYFASTTREE {
    tag "$alignment"
    label 'process_iqtree'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/veryfasttree_3.2.1:0.0.1':
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/veryfasttree_3.2.1:0.0.1':
        '' }"

    publishDir(
        path: "${params.outdir}/veryfasttree_gene_trees",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(alignment)
    val model // not used

    output:
    path("*.treefile")  , emit: phylogeny
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    """
    og=\${alignment%%_*}
 
    # Infer a... Very Fast Tree!
    VeryFastTree \\
        -threads ${task.cpus} \\
        $alignment > \${og}_vft.treefile
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VeryFastTree: $(VeryFastTree | head -n1 | cut -f2 -d" ")
    END_VERSIONS
    """
}
