process TRANSDECODER {
    tag "${meta.id}"
    label 'process_medium'

    container 'arcadiascience/preprocess_proteomes:1.1.0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_translated.fasta"), emit: translated
    path "versions.yml",                                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Decompress if gzipped
    input_file="${fasta}"
    if [[ "${fasta}" == *.gz ]]; then
        gunzip -c "${fasta}" > input.fasta
        input_file="input.fasta"
    fi

    # Detect if input is already protein sequences by sampling non-header lines.
    # Nucleotide sequences should be almost entirely ACGTN; protein sequences
    # contain amino acid letters (D, E, F, H, I, K, L, M, P, Q, R, W, Y) that
    # are absent or very rare in nucleotide data.
    aa_chars=\$(grep -v '^>' "\${input_file}" | head -1000 | tr -d '\\n ' | grep -o '[DEFHIKLMPQRWYdefhiklmpqrwy]' | wc -c || echo 0)
    total_chars=\$(grep -v '^>' "\${input_file}" | head -1000 | tr -d '\\n ' | wc -c || echo 0)

    if [ "\${total_chars}" -gt 0 ]; then
        aa_frac=\$(awk "BEGIN {printf \\"%.2f\\", \${aa_chars}/\${total_chars}}")
    else
        aa_frac="0.00"
    fi

    # If >20% of characters are amino-acid-only letters, this is protein — skip TransDecoder
    is_protein=\$(awk "BEGIN {print (\${aa_frac} > 0.20) ? 1 : 0}")

    if [ "\${is_protein}" -eq 1 ]; then
        echo "WARNING: Input for ${meta.id} appears to be protein sequences (aa fraction: \${aa_frac}). Skipping TransDecoder." >&2
        cp "\${input_file}" "${meta.id}_translated.fasta"
    else
        TransDecoder.LongOrfs -t "\${input_file}"
        TransDecoder.Predict -t "\${input_file}" --single_best_only
        cp "\${input_file}.transdecoder.pep" "${meta.id}_translated.fasta"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transdecoder: \$(TransDecoder.LongOrfs --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo "5.7.1")
    END_VERSIONS
    """
}
