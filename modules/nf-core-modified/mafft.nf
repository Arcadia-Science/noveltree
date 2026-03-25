process MAFFT {
    // Modified from nf-core to:
    // 1) Update input/output to match tuple formatting used throughout workflow
    // 2) use MAFFTs built in automatic thread scaling to improve memory efficiency
    // 3) create protein/species map file within module
    tag "${meta.og}"

    cpus { Math.min( 12 * task.attempt, params.max_cpus as int ) }
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def library_gb = n * n * L * 16L / (1024L * 1024L * 1024L)
        def dp_gb = L * L * 24L / (1024L * 1024L * 1024L)
        def estimated_gb = Math.max(8L, (long)(library_gb + dp_gb) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 128L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    conda (params.enable_conda ? "bioconda::mafft=7.490" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.490--h779adbc_0':
        'quay.io/biocontainers/mafft:7.490--h779adbc_0' }"

    storeDir "${params.outdir}/alignments/original"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}_mafft.fa")          , emit: msas
    tuple val(meta), path("species_protein_maps/${fasta.baseName}_map.link"), emit: map_link, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def aln_trimmer = params.msa_trimmer
    """
    # Be sure to remove any non-standard amino acid codes in the input sequences, as this
    # can cause errors downstream and in parsing.
    sed -E -i '/>/!s/U/X/g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O/X/g' ${fasta} # pyrrolysine

    prefix=\$(basename "${fasta}" .fa)
    mafft \\
        --thread ${task.cpus} \\
        ${args} ${fasta} > \${prefix}_mafft.fa

    # Create protein-species map files if we are not doing any alignment cleaning
    if [ $aln_trimmer == "none" ]; then
        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        mkdir species_protein_maps
        grep ">" \${prefix}_mafft.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > species_protein_maps/\${prefix}_map.link
        rm prot && rm spp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
