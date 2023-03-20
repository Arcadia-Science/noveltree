process FASTTREE {
    tag "$alignment"
    label 'process_iqtree'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/magus_0.1.0:0.0.1':
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/magus_0.1.0:0.0.1':
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
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    og=\$(echo $alignment | cut -f1 -d "_")
    
    # Make sure the number of threads are being specified properly 
    export OMP_NUM_THREADS=${task.cpus}
    
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
