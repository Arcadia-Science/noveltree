process PREPROCESS_PROTEOME {
    tag "${meta.id}"
    label 'process_medium'

    container 'arcadiascience/preprocess_proteomes:1.1.0'

    storeDir "${params.outdir}/preprocessing/preprocessed_proteomes"

    input:
    tuple val(meta), path(fasta)
    val min_protein_length

    output:
    tuple val(meta), path("${meta.id}_preprocessed.fasta"), emit: preprocessed

    when:
    task.ext.when == null || task.ext.when

    script:
    def skip_cdhit = meta.reference == 'yes'
    def cdhit_threshold = meta.transdecoder == 'yes' ? '0.97' : '1.00'
    """
    # Decompress if gzipped
    input_file="${fasta}"
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > input.fasta
        input_file="input.fasta"
    fi

    # Remove trailing stop codons (*), replace internal stops with X
    # Handle rare amino acids: U→C (selenocysteine), J/B/Z→X (ambiguous)
    sed '/^>/!{
        s/\\*\$//
        s/\\*/X/g
        s/U/C/g
        s/J/X/g
        s/B/X/g
        s/Z/X/g
    }' "\${input_file}" > cleaned.fasta

    # Filter by minimum length
    seqkit seq --min-len ${min_protein_length} cleaned.fasta > length_filtered.fasta

    # CD-HIT clustering (conditional on meta flags)
    if [ "${skip_cdhit}" == "true" ]; then
        # Reference proteomes: skip CD-HIT (already one protein per gene)
        cp length_filtered.fasta "${meta.id}_preprocessed.fasta"
    else
        cd-hit -i length_filtered.fasta -o "${meta.id}_preprocessed.fasta" \
            -c ${cdhit_threshold} -n 5 -M 0 -T ${task.cpus}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/seqkit v//')
        cd-hit: \$(cd-hit -h 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo "4.8.1")
    END_VERSIONS
    """
}
