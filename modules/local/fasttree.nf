process FASTTREE {
    tag "$meta.og"

    cpus {
        def n = (meta?.n_seq ?: 0) as long
        n >= (params.fasttree_large_family_min_sequences as long) \
            ? Math.min(
                (params.fasttree_large_family_cpus as int) +
                    ((task.attempt - 1) * (params.fasttree_large_family_retry_cpu_step as int)),
                params.max_cpus as int
            ) \
            : Math.min(16 * task.attempt, params.max_cpus as int)
    }
    time {
        def n = (meta?.n_seq ?: 0) as long
        def maxTime = params.max_time as nextflow.util.Duration
        def requested = n >= (params.fasttree_large_family_min_sequences as long) \
            ? (task.attempt == 1 \
                ? params.fasttree_large_family_time as nextflow.util.Duration \
                : params.fasttree_large_family_retry_time as nextflow.util.Duration) \
            : 3.h * Math.pow(3, task.attempt - 1)
        requested.compareTo(maxTime) > 0 ? maxTime : requested
    }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def profile_bytes = n * L * 165L * 3L  // 3x multiplier for SPR/NNI search overhead
        def nj_bytes = (long)(16.0 * Math.pow(n, 1.5))
        def estimated_gb = Math.max(8L, (long)((profile_bytes + nj_bytes) / (1024L * 1024L * 1024L)) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 128L)
        def multiplier = n >= (params.fasttree_large_family_min_sequences as long) \
            ? (params.fasttree_large_family_memory_multiplier as int) + task.attempt - 1 \
            : task.attempt
        def requested = capped_gb.GB * multiplier
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/fasttree_2.1.11:1.0.0'

    storeDir "${params.outdir}/gene_family_trees/original"

    input:
    tuple val(meta), file(alignment)
    val model // not used

    output:
    tuple val(meta), path("${alignment.baseName}_ft.newick") , emit: phylogeny

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        FastTree: \$(FastTreeDblMP 2>&1 | head -n1 | cut -d" " -f5)
    END_VERSIONS
    """
}
