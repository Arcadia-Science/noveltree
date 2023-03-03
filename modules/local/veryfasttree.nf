process VERYFASTTREE {
    tag "$alignment"
    label 'process_iqtree'

    container "${ workflow.containerEngine == 'docker' ? 'quay.io/biocontainers/veryfasttree:3.2.1--h9f5acd7_0':
        '' }"

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
    
    # Hacky fix to prevent segfault of VeryFastTree when running on small datasets
    nseqs=\$(grep ">" $alignment | wc -l)
    if [[ \${nseqs} -gt 14 ]]
    then
        nthreads=${task.cpus}
    else
        nthreads=1
    fi
    
    # Infer a... Very Fast Tree!
    VeryFastTree \\
        -threads \${nthreads} \\
        $alignment > \${og}_vft.treefile
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VeryFastTree: \$(VeryFastTree | head -n1 | cut -f2 -d" ")
    END_VERSIONS
    """
}
