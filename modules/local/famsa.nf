process FAMSA {
    tag "$meta.og"

    cpus { Math.min( 36 * task.attempt, params.max_cpus as int ) }
    time { 12.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 5000) as long
        def L = (meta?.max_len ?: 500) as long
        def L_aln = L * 3L
        def nCpus = task.cpus as long
        // Progressive alignment: profiles + parallel DP matrices (one per concurrent merge)
        def profile_bytes = n * L_aln * 96L          // 24 symbols × 4 bytes
        def dp_bytes = L_aln * L_aln * 16L
        def parallel_merges = (long) Math.min(nCpus, Math.max(1L, (long) Math.sqrt(n as double)))
        def estimated_gb = Math.max(32L, (long)((profile_bytes + dp_bytes * parallel_merges) / (1024L * 1024L * 1024L)) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 256L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/famsa_2.0.0:1.0.0'

    storeDir "${params.outdir}/alignments/original"

    input:
    tuple val(meta), path(fasta), path(family_map)

    output:
    tuple val(meta), path("${fasta.baseName}_famsa.fa")          , emit: msas
    tuple val(meta), path("species_protein_maps/${fasta.baseName}_map.link"), emit: map_link, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = fasta.baseName
    """
    # Run FAMSA alignment
    famsa \\
        -t ${task.cpus} \\
        ${args} \\
        ${fasta} \\
        ${prefix}_famsa.fa

    mkdir species_protein_maps
    cp ${family_map} species_protein_maps/${prefix}_map.link

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famsa: \$(famsa 2>&1 | head -n 1 | sed 's/.*FAMSA //' | sed 's/ .*//')
    END_VERSIONS
    """
}
