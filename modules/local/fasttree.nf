process FASTTREE {
    tag "$meta.og"

    cpus { Math.min( 16 * task.attempt, params.max_cpus as int ) }
    time { 3.h * Math.pow(3, task.attempt - 1) }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def profile_bytes = n * L * 165L * 3L  // 3x multiplier for SPR/NNI search overhead
        def nj_bytes = (long)(16.0 * Math.pow(n, 1.5))
        def estimated_gb = Math.max(8L, (long)((profile_bytes + nj_bytes) / (1024L * 1024L * 1024L)) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 128L)
        def requested = capped_gb.GB * task.attempt
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
