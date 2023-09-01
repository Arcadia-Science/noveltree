process CIALIGN {
    tag "$fasta"
    label 'process_lowcpu'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/cialign_1.1.0:1.0.0' :
        '' }"

    publishDir(
        path: "${params.outdir}/cialign_cleaned_msas",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)              // Filepaths to the MSAs

    output:
    tuple val(meta), path("**_cialign.fa") , emit: cleaned_msas, optional: true
    tuple val(meta), path("**_map.link")   , emit: map_link, optional: true
    path "*"                               , emit: results
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ''
    def remove_short = params.min_ungapped_length ? "--remove_short --remove_min_length=${params.min_ungapped_length}" : ''
    """
    # Get the name of the orthogroup we are processing
    prefix=\$(echo ${fasta} | cut -f1 -d "_")

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

    # Now, create a protein-species map-file, assuming that trimming didn't lead
    # to the focal MSA being comprised of < 4 sequences.
    n_remain=\$(grep ">" \${prefix}_cialign.fa | wc -l)
    if [ \$n_remain -lt 4 ]; then
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
