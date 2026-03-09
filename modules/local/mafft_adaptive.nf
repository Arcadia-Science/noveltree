process MAFFT_ADAPTIVE {
    tag "$meta.og"

    cpus { 12 * task.attempt }
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def library_gb = n * n * L * 16L / (1024L * 1024L * 1024L)
        def dp_gb = L * L * 24L / (1024L * 1024L * 1024L)
        def estimated_gb = Math.max(4L, (long)(library_gb + dp_gb) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 48L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    conda (params.enable_conda ? "bioconda::mafft=7.490" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.490--h779adbc_0':
        'quay.io/biocontainers/mafft:7.490--h779adbc_0' }"

    publishDir(
        path: "${params.outdir}/alignments/original",
        mode: params.publish_dir_mode,
        pattern: "*_{einsi,linsi}.fa",
    )
    publishDir(
        path: "${params.outdir}/alignments/species_protein_maps",
        mode: params.publish_dir_mode,
        pattern: "species_protein_maps/*",
        saveAs: { fn -> fn.split('/')[-1] },
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_{einsi,linsi}.fa")    , emit: msas
    tuple val(meta), path("species_protein_maps/*_map.link"), emit: map_link, optional: true
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def aln_trimmer = params.msa_trimmer
    def mafft_mode  = params.mafft_mode ?: 'einsi'
    def mafft_args  = mafft_mode == 'linsi' ?
        '--localpair --maxiterate 1000 --anysymbol' :
        '--genafpair --maxiterate 1000 --ep 0 --anysymbol'
    def prefix = fasta.baseName
    """
    # Remove non-standard amino acid codes
    sed -E -i '/>/!s/U/X/g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O/X/g' ${fasta} # pyrrolysine

    mafft \\
        --thread ${task.cpus} \\
        ${mafft_args} \\
        ${fasta} > ${prefix}_${mafft_mode}.fa

    # Create protein-species map files if we are not doing any alignment cleaning
    if [ "${aln_trimmer}" == "none" ]; then
        mkdir species_protein_maps
        grep ">" ${prefix}_${mafft_mode}.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > species_protein_maps/${prefix}_map.link
        rm prot && rm spp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
    END_VERSIONS
    """
}
