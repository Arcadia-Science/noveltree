process MAFFT_TIER1 {
    tag "$meta.og"

    cpus { Math.min( 12 * task.attempt, params.max_cpus as int ) }
    time { 12.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def library_gb = n * n * L * 16L / (1024L * 1024L * 1024L)
        def dp_gb = L * L * 24L / (1024L * 1024L * 1024L)
        def estimated_gb = Math.max(32L, (long)(library_gb + dp_gb) + 4L)
        def capped_gb = (int) Math.min(estimated_gb, 256L)
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
    tuple val(meta), path(fasta), path(family_map)

    output:
    tuple val(meta), path("${fasta.baseName}_${(meta.n_seq as int) <= (params.align_tier1_max as int) ? 'einsi' : 'linsi'}.fa"), emit: msas
    tuple val(meta), path("species_protein_maps/${fasta.baseName}_map.link"), emit: map_link, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def mafft_mode  = params.mafft_mode ?: 'einsi'
    def mafft_args  = mafft_mode == 'linsi' ?
        '--localpair --maxiterate 1000 --anysymbol' :
        '--genafpair --maxiterate 1000 --ep 0 --anysymbol'
    def prefix = fasta.baseName
    """
    mafft \\
        --thread ${task.cpus} \\
        ${mafft_args} \\
        ${fasta} > ${prefix}_${mafft_mode}.fa

    mkdir species_protein_maps
    cp ${family_map} species_protein_maps/${prefix}_map.link

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
