process CLIPKIT {
    tag "$fasta"
    label 'process_low_cpu'

    container 'arcadiascience/clipkit_2.1.1-seqmagick_0.8.4:1.0.0'

    publishDir(
        path: "${params.outdir}/alignments/trimmed",
        mode: params.publish_dir_mode,
        pattern: "*_clipkit.fa",
    )
    publishDir(
        path: "${params.outdir}/alignments/species_protein_maps",
        mode: params.publish_dir_mode,
        pattern: "species_protein_maps/*",
        saveAs: { fn -> fn.split('/')[-1] },
    )

    input:
    tuple val(meta), path(fasta)              // Filepaths to the MSAs

    output:
    tuple val(meta), path("*_clipkit.fa")  , emit: cleaned_msas, optional: true
    tuple val(meta), path("species_protein_maps/*_map.link"), emit: map_link, optional: true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    def min_ungapped_length = params.min_ungapped_length
    def min_seq = params.min_num_seq_per_og
    def min_spp = params.min_num_spp_per_og
    """
    # Get the alignment prefix (strip .fa extension, preserving aligner provenance)
    prefix=\$(basename "$fasta" .fa)

    # Trim the MSAs for each orthogroup containing at least 4 species.
    clipkit ${fasta} -o \${prefix}_tmp.fa $args

    # Remove sequences with a minimum non-gapped length less than the specified length.
    seqmagick convert \\
        --min-ungapped-length $min_ungapped_length \\
        \${prefix}_tmp.fa \\
        \${prefix}_clipkit.fa

    # Verify the trimmed alignment still meets minimum sequence/species thresholds.
    # Trimming can remove sequences, potentially dropping an OG below the filters
    # that were applied to the raw FASTA files.
    n_seq=\$(grep -c ">" \${prefix}_clipkit.fa || true)
    n_spp=\$(grep ">" \${prefix}_clipkit.fa | sed "s/>//" | sed "s/_[^_]*\$//" | sort -u | wc -l | tr -d ' ')
    if [ "\$n_seq" -lt "$min_seq" ] || [ "\$n_spp" -lt "$min_spp" ]; then
        rm \${prefix}_clipkit.fa
    else
        # Now, create a protein-species map-file:
        # Pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        mkdir species_protein_maps
        grep ">" \${prefix}_clipkit.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > species_protein_maps/\${prefix}_map.link
        rm prot && rm spp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clipkit: \$( clipkit --version | sed "s/clipkit //g" )
        seqmagick: \$( seqmagick --version | cut -f2 -d" " )
    END_VERSIONS
    """
}
