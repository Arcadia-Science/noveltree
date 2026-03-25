process PARSE_PHYLOHOGS {
    tag "$meta.og"

    cpus 1
    time { 4.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 100) as long
        // Tree traversal is O(n); pair output is O(n²) but streamed to disk
        def base_gb = Math.max(4L, (long)(n / 200L) + 4L)
        def capped_gb = (int) Math.min(base_gb, 64L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/phylo_profiles:1.0.0'

    storeDir "${params.outdir}/orthology"

    input:
    tuple val(meta), path(nhx_file)
    path species_tree

    output:
    tuple val(meta), path("${meta.og}/${meta.og}_orthologs.tsv") , emit: orthologs
    tuple val(meta), path("${meta.og}/${meta.og}_paralogs.tsv")  , emit: paralogs
    path "${meta.og}/${meta.og}_hog_membership.tsv"    , emit: hog_membership
    path "spp_tree_node_lookup.tsv"                   , emit: node_lookup

    when:
    task.ext.when == null || task.ext.when

    script:
    def og = "${meta.og}"
    """
    mkdir -p ${og}
    extract_relationships_from_nhx.py ${nhx_file} ${species_tree} ${og} ${og}/${og}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | cut -d' ' -f2 )
    END_VERSIONS
    """
}
