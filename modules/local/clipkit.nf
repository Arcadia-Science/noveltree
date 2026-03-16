process CLIPKIT {
    tag "${meta.og}"

    cpus 1
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def estimated_gb = Math.max(2L, (long)(n * L * 20L / (1024L * 1024L * 1024L)) + 1L)
        def capped_gb = (int) Math.min(estimated_gb, 32L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/clipkit_2.1.1-seqmagick_0.8.4:1.0.0'

    storeDir "${params.outdir}/alignments/trimmed"

    input:
    tuple val(meta), path(fasta)              // Filepaths to the MSAs

    output:
    tuple val(meta), path("${fasta.baseName}_clipkit.fa")  , emit: cleaned_msas
    tuple val(meta), path("species_protein_maps/${fasta.baseName}_map.link"), emit: map_link

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

    # Check if the trimmed alignment has fewer columns than min_ungapped_length.
    # If so, no sequence can pass the filter — skip seqmagick entirely.
    n_cols=\$(awk '!/^>/{print length; exit}' \${prefix}_tmp.fa)
    if [ "\$n_cols" -lt "$min_ungapped_length" ]; then
        # Alignment too short after trimming — will produce empty outputs below
        rm \${prefix}_tmp.fa
    else
        # Remove sequences with a minimum non-gapped length less than the specified length.
        seqmagick convert \\
            --min-ungapped-length $min_ungapped_length \\
            \${prefix}_tmp.fa \\
            \${prefix}_clipkit.fa
    fi

    # Verify the trimmed alignment still meets minimum sequence/species thresholds.
    # Trimming can remove sequences, potentially dropping an OG below the filters
    # that were applied to the raw FASTA files.
    if [ -f \${prefix}_clipkit.fa ]; then
        n_seq=\$(grep -c ">" \${prefix}_clipkit.fa || true)
        n_spp=\$(grep ">" \${prefix}_clipkit.fa | sed "s/>//" | sed "s/_[^_]*\$//" | sort -u | wc -l | tr -d ' ')
    else
        n_seq=0
        n_spp=0
    fi
    mkdir -p species_protein_maps
    if [ "\$n_seq" -lt "$min_seq" ] || [ "\$n_spp" -lt "$min_spp" ]; then
        # QC failed — produce empty outputs so storeDir can distinguish
        # "task ran, alignment discarded" from "task never ran"
        rm -f \${prefix}_clipkit.fa
        touch \${prefix}_clipkit.fa
        touch species_protein_maps/\${prefix}_map.link
    else
        # Create a protein-species map-file:
        # Pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
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
