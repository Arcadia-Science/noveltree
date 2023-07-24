process CLIPKIT {
    tag "$fasta"
    label 'process_lowcpu'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/clipkit_1.3.0-seqmagick_0.8.4:0.0.1' :
        '' }"

    publishDir(
        path: "${params.outdir}/clipkit_cleaned_msas",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)              // Filepaths to the MSAs

    output:
    tuple val(meta), path("**_clipkit.fa") , emit: cleaned_msas, optional: true
    tuple val(meta), path("**_map.link")   , emit: map_link, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    def min_ungapped_length = params.min_ungapped_length
    """
    # Get the name of the orthogroup we are processing
    prefix=\$(echo $fasta | cut -f1 -d "_")

    # Trim the MSAs for each orthogroup containing at least 4 species.
    clipkit ${fasta} -o \${prefix}_tmp.fa $args

    # Remove sequences with a minimum non-gapped length less than the specified length.
    seqmagick convert \\
        --min-ungapped-length $min_ungapped_length \\
        \${prefix}_tmp.fa \\
        \${prefix}_clipkit.fa

    # Now, create a protein-species map-file: 
    # Pull out the sequences, and split into a TreeRecs format mapping
    # file, where each protein in the tree is a new line, listing species
    # and then the protein
    mkdir species_protein_maps
    grep ">" \${prefix}_clipkit.fa | sed "s/>//g"  | sed "s/.*://g" > prot
    sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
    paste prot spp > species_protein_maps/\${prefix}_map.link
    rm prot && rm spp
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$( clipkit --version | sed "s/clipkit //g" )
        seqmagick: \$( seqmagick --version | cut -f2 -d" " )
    END_VERSIONS
    """
}
