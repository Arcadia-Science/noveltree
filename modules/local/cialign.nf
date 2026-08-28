process CIALIGN {
    tag "$fasta"
    label 'process_low_cpu'

    container 'arcadiascience/cialign_1.1.0:1.0.0'

    storeDir "${params.outdir}/alignments/trimmed"

    input:
    tuple val(meta), path(fasta), path(family_map)

    output:
    tuple val(meta), path("${fasta.baseName}_cialign.fa")  , emit: cleaned_msas
    tuple val(meta), path("species_protein_maps/${fasta.baseName}_map.link"), emit: map_link
    path "removed_sites/*"                 , emit: removed_sites
    path "log_files/*"                     , emit: log_files

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
    if [ "\$n_seq" -gt 0 ]; then
        subset_gene_species_map.py \
            --mapping ${family_map} \
            --fasta \${prefix}_cialign.fa \
            --output \${prefix}_map.link
        n_spp=\$(cut -f2 \${prefix}_map.link | sort -u | wc -l | tr -d ' ')
    else
        n_spp=0
    fi
    mkdir -p species_protein_maps
    if [ "\$n_seq" -lt "$min_seq" ] || [ "\$n_spp" -lt "$min_spp" ]; then
        # QC failed — produce empty outputs so storeDir can distinguish
        # "task ran, alignment discarded" from "task never ran"
        rm \${prefix}_cialign.fa
        touch \${prefix}_cialign.fa
        touch species_protein_maps/\${prefix}_map.link
    else
        mv \${prefix}_map.link species_protein_maps/\${prefix}_map.link
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CIAlign: \$( CIAlign --version )
    END_VERSIONS
    """
}
