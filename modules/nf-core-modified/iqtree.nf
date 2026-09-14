process IQTREE {
    // Modified from nf-core to:
    // 1) remove constant sites specification (not applicable to workflow)
    // 2) parameterize tree-model specification
    // 3) correctly handle "task.memory" specification for memory handling by iqtree
    // 4) update the Docker container to use iqtree v2.2.0.5
    // 6) input/output files in appropriate tuple format
    tag "${meta.og}"

    cpus { Math.min( 12, params.max_cpus as int ) }
    time {
        params.iqtree_fasttree_fallback \
            ? (params.iqtree_timeout as nextflow.util.Duration) * (Math.pow(3, task.attempt - 1) as int) +
                (params.iqtree_timeout_grace as nextflow.util.Duration) \
            : 3.h * Math.pow(3, task.attempt - 1)
    }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def clv_bytes = n * L * 640L
        def estimated_gb = Math.max(16L, (long)(clv_bytes * 1.5 / (1024L * 1024L * 1024L)) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 128L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/iqtree_2.2.0.5:1.0.0'

    storeDir "${params.outdir}/gene_family_trees/original"

    input:
    tuple val(meta), file(alignment)
    val model

    output:
    tuple val(meta), path("${alignment.baseName}_iqt.{newick,fallback.json}", arity: '1'), emit: outcome
    tuple val(meta), path("${alignment.baseName}_iqt.newick"), emit: phylogeny, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def memory   = task.memory.toString().replaceAll(' ', '')
    def prefix   = alignment.baseName
    def runner = ''
    if (params.iqtree_fasttree_fallback) {
        def budget = (params.iqtree_timeout as nextflow.util.Duration) * (Math.pow(3, task.attempt - 1) as int)
        def grace = params.iqtree_timeout_grace as nextflow.util.Duration
        if (budget.toMillis() < 1000 || grace.toMillis() < 2000 || task.time < budget + grace) {
            throw new IllegalArgumentException('IQ-TREE needs a positive computation budget and at least 2 seconds of grace; task.time must include both')
        }
        def killAfter = Math.min(30L, Math.max(1L, grace.toSeconds().intdiv(2)))
        runner = "bash run_iqtree_with_fallback.sh '${meta.og}' ${budget.toSeconds()} ${killAfter} '${alignment}.treefile' '${prefix}_iqt.newick' '${prefix}_iqt.fallback.json' --"
    }
    """
    memory=\$(echo ${task.memory} | sed "s/.G/G/g")

    # Infer the phylogeny
    test -s "$alignment"
    ${runner} iqtree2 \\
        -s $alignment \\
        -nt AUTO \\
        -ntmax ${task.cpus} \\
        -mem \$memory \\
        -m $model \\
        $args

    # Rename to standardized output format
    ${params.iqtree_fasttree_fallback ? '' : "mv '${alignment}.treefile' '${prefix}_iqt.newick'"}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(iqtree2 -version 2>&1 | head -n1)
    END_VERSIONS
    """
}
