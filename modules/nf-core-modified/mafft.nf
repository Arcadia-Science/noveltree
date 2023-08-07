process MAFFT {
    // Modified from nf-core to:
    // 1) Update input/output to match tuple formatting used throughout workflow
    // 2) use MAFFTs built in automatic thread scaling to improve memory efficiency
    // 3) create protein/species map file within module
    tag "$fasta"
    label 'process_mafft'

    conda (params.enable_conda ? "bioconda::mafft=7.490" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.490--h779adbc_0':
        'quay.io/biocontainers/mafft:7.490--h779adbc_0' }"

    publishDir(
        path: "${params.outdir}/mafft_alignments",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("**_mafft.fa")         , emit: msas
    tuple val(meta), path("**_map.link")         , emit: map_link, optional: true
    path("*")                                    , emit: results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def aln_trimmer = params.msa_trimmer
    """
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
