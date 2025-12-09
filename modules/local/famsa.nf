process FAMSA {
    tag "$meta.og"
    label 'process_high'

    container 'arcadiascience/famsa_2.0.0:1.0.0'

    publishDir(
        path: "${params.outdir}/famsa_alignments",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("**_famsa.fa")         , emit: msas
    tuple val(meta), path("**_map.link")         , emit: map_link, optional: true
    path("*")                                    , emit: results
    path "versions.yml"                          , emit: versions

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
