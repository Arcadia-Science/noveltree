process IQTREE {
    // Modified from nf-core to:
    // 1) remove constant sites specification (not applicable to workflow)
    // 2) parameterize tree-model specification
    // 3) correctly handle "task.memory" specification for memory handling by iqtree
    // 4) update the Docker container to use iqtree v2.2.0.5
    tag "$alignment"
    label 'process_iqtree'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/iqtree_2.2.0.5:0.0.1':
        '' }"

    input:
    path(alignment)
    val model

    output:
    path("*.treefile")  , emit: phylogeny
    path(alignment)     , emit: msa
    path "*.log"        , emit: iqtree_log
    path "versions.yml" , emit: versions

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
