process IQTREE_PMSF {
    // Modified from nf-core to:
    // 1) remove constant sites specification (not applicable to workflow)
    // 2) parameterize tree-model specification
    // 3) optionally infer trees using PMSF approximation which requires initial tree inference
    // 4) correctly handle "task.memory" specification for memory handling by iqtree
    // 5) update the Docker container to use iqtree v2.2.0.5
    tag "$alignment"
    label 'process_iqtree'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/iqtree_2.2.0.5:0.0.1':
        '' }"
    
    publishDir(
        path: "${params.outdir}/iqtree_pmsf_gene_trees",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(alignment)
    path(guide_tree)
    path(guide_tree_log)
    val pmsf_model

    output:
    path("*pmsf.treefile") , emit: phylogeny
    path "*.log"           , emit: iqtree_log
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def memory = task.memory.toString().replaceAll(' ', '')

    """
    memory=\$(echo ${task.memory} | sed "s/.G/G/g")

    # Identify the best number of threads
    nt=\$(grep "BEST NUMBER" $guide_tree_log | sed "s/.*: //g")

    # Rename things and clean up - iqtree will be unhappy otherwise
    mv $guide_tree guidetree.treefile
    rm $guide_tree_log

    iqtree2 \\
        -s $alignment \\
        -nt \$nt \\
        -mem \$memory \\
        -m $pmsf_model \\
        -f $guide_tree \\
        -ft $guide_tree \\
        $args
    
    # Clean up
    rm guidetree.treefile
    
    # Rename the tree and log file to something more informative
    mv ${alignment}.treefile ${alignment}.pmsf.treefile
    mv ${alignment}.log ${alignment}.pmsf.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
