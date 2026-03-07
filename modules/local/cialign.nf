process CIALIGN {
    tag "$fasta"
    label 'process_low_cpu'

    container 'arcadiascience/cialign_1.1.0:1.0.0'

    publishDir(
        path: "${params.outdir}/alignments/trimmed",
        mode: params.publish_dir_mode,
        pattern: "*_cialign.fa",
    )
    publishDir(
        path: "${params.outdir}/alignments/species_protein_maps",
        mode: params.publish_dir_mode,
        pattern: "species_protein_maps/*",
        saveAs: { fn -> fn.split('/')[-1] },
    )
    publishDir(
        path: "${params.outdir}/alignments/trimmed",
        mode: params.publish_dir_mode,
        pattern: "{removed_sites,log_files}/*",
        saveAs: { fn -> fn },
    )

    input:
    tuple val(meta), path(fasta)              // Filepaths to the MSAs

    output:
    tuple val(meta), path("*_cialign.fa")  , emit: cleaned_msas, optional: true
    tuple val(meta), path("species_protein_maps/*_map.link"), emit: map_link, optional: true
    path "removed_sites/*"                 , emit: removed_sites
    path "log_files/*"                     , emit: log_files
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def remove_short = params.min_ungapped_length ? "--remove_short --remove_min_length=${params.min_ungapped_length}" : ''
    def min_seq = params.min_num_seq_per_og
    def min_spp = params.min_num_spp_per_og
    """
    # Get the alignment prefix (strip .fa extension, preserving aligner provenance)
    prefix=\$(basename "${fasta}" .fa)

    # Clean up the MSAs for each orthogroup containing at least 4 species.
    CIAlign \
        --infile ${fasta} \
        --outfile_stem="\${prefix}" \
        ${remove_short} \
        $args

    # Rename output so it is clear we trimmed with CIAlign
    mv \${prefix}_cleaned.fasta \${prefix}_cialign.fa

    # And move the "removed.txt" files indicating which sites were removed
    # from each MSA to a separate directory
    mkdir -p removed_sites
    mv *removed.txt removed_sites

    # And do the same for the log files
    mkdir log_files
    mv *log.txt log_files

    # Verify the trimmed alignment still meets minimum sequence/species thresholds.
    n_seq=\$(grep -c ">" \${prefix}_cialign.fa || true)
    n_spp=\$(grep ">" \${prefix}_cialign.fa | sed "s/>//" | sed "s/_[^_]*\$//" | sort -u | wc -l | tr -d ' ')
    if [ "\$n_seq" -lt "$min_seq" ] || [ "\$n_spp" -lt "$min_spp" ]; then
        rm \${prefix}_cialign.fa
    else
        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        mkdir species_protein_maps
        grep ">" \${prefix}_cialign.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > species_protein_maps/\${prefix}_map.link
        rm prot && rm spp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CIAlign: \$( CIAlign --version )
    END_VERSIONS
    """
}
