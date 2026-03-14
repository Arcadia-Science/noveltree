process FAMSA {
    tag "$meta.og"

    cpus { 36 * task.attempt }
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 5000) as long
        def L = (meta?.max_len ?: 500) as long
        def L_aln = L * 3L
        def nCpus = task.cpus as long
        // Progressive alignment: profiles + parallel DP matrices (one per concurrent merge)
        def profile_bytes = n * L_aln * 96L          // 24 symbols × 4 bytes
        def dp_bytes = L_aln * L_aln * 16L
        def parallel_merges = (long) Math.min(nCpus, Math.max(1L, (long) Math.sqrt(n as double)))
        def estimated_gb = Math.max(4L, (long)((profile_bytes + dp_bytes * parallel_merges) / (1024L * 1024L * 1024L)) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 128L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/famsa_2.0.0:1.0.0'

    storeDir "${params.outdir}/alignments/original"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}_famsa.fa")          , emit: msas
    tuple val(meta), path("species_protein_maps/${fasta.baseName}_map.link"), emit: map_link, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def aln_trimmer = params.msa_trimmer
    def prefix = fasta.baseName
    """
    # Be sure to remove any non-standard amino acid codes in the input sequences, as this
    # can cause errors downstream and in parsing.
    sed -E -i '/>/!s/U/X/g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O/X/g' ${fasta} # pyrrolysine

    # Run FAMSA alignment
    famsa \\
        -t ${task.cpus} \\
        ${args} \\
        ${fasta} \\
        ${prefix}_famsa.fa

    # Create protein-species map files if we are not doing any alignment cleaning
    if [ "${aln_trimmer}" == "none" ]; then
        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        mkdir species_protein_maps
        grep ">" ${prefix}_famsa.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > species_protein_maps/${prefix}_map.link
        rm prot && rm spp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        famsa: \$(famsa 2>&1 | head -n 1 | sed 's/.*FAMSA //' | sed 's/ .*//')
    END_VERSIONS
    """
}
