process IQTREE {
    // Modified from nf-core to:
    // 1) remove constant sites specification (not applicable to workflow)
    // 2) parameterize tree-model specification
    // 3) correctly handle "task.memory" specification for memory handling by iqtree
    // 4) update the Docker container to use iqtree v2.2.0.5
    // 6) input/output files in appropriate tuple format
    tag "$alignment"
    label 'process_iqtree'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/iqtree_2.2.0.5:1.0.0':
        '' }"

    publishDir(
        path: "${params.outdir}/iqtree_gene_trees",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), file(alignment)
    val model

    output:
    tuple val(meta), path("*.treefile") , emit: phylogeny
    tuple val(meta), path("*.log")      , emit: iqtree_log
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def memory = task.memory.toString().replaceAll(' ', '')

    """
    memory=\$(echo ${task.memory} | sed "s/.G/G/g")

    # Check if this is a resumed run:
    # error trying to resume if not.)
    # If the checkpoint file indicates the run finished, go ahead and
    # skip the analyses, otherwise run iqtree as normal.

    # Infer the phylogeny
    iqtree2 \\
        -s $alignment \\
        -nt AUTO \\
        -ntmax ${task.cpus} \\
        -mem \$memory \\
        -m $model \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """
}
