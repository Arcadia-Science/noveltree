process DOWNLOAD_INPUT {
    tag "${meta.id}"
    label 'process_single'

    container 'arcadiascience/preprocess_proteomes:1.1.0'

    input:
    tuple val(meta), val(source)

    output:
    tuple val(meta), path("${meta.id}_downloaded.*"), emit: downloaded
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (meta.source_type == 'ncbi_refseq') {
        """
        datasets download genome accession ${source} --include protein --filename data.zip
        unzip -o data.zip -d ncbi_data

        PROTEIN_FILE=\$(find ncbi_data -name 'protein.faa' | head -1)
        if [ -z "\${PROTEIN_FILE}" ]; then
            echo "ERROR: No protein.faa found for accession ${source}. This genome may lack protein annotations." >&2
            exit 1
        fi

        mv "\${PROTEIN_FILE}" "${meta.id}_downloaded.fasta"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            datasets: \$(datasets --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo "unknown")
        END_VERSIONS
        """
    } else if (meta.source_type == 'uniprot') {
        """
        curl -sL "https://rest.uniprot.org/uniprotkb/stream?query=(proteome:${source})&format=fasta" \
            -o "${meta.id}_downloaded.fasta"

        if [ ! -s "${meta.id}_downloaded.fasta" ]; then
            echo "ERROR: Empty result for UniProt proteome ${source}. Check that the proteome ID is valid." >&2
            exit 1
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            curl: \$(curl --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' || echo "unknown")
        END_VERSIONS
        """
    } else if (meta.source_type == 'ncbi_tsa') {
        // TSA accession — construct WGS/TSA FTP URL and download nucleotide FASTA
        def chars = source.replaceAll(/\d+$/, '')  // e.g. GKVO01000000 -> GKVO
        def prefix = source[0..5]  // e.g. GKVO01
        def d1 = source[0..1]      // e.g. GK
        def d2 = source[2..3]      // e.g. VO
        """
        wget -q -O "${meta.id}_downloaded.fna.gz" \
            "https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/${d1}/${d2}/${prefix}/${prefix}.1.fsa_nt.gz"

        if [ ! -s "${meta.id}_downloaded.fna.gz" ]; then
            echo "ERROR: Download failed for TSA accession ${source}. Check that the accession is valid." >&2
            exit 1
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            wget: \$(wget --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+' || echo "1.21")
        END_VERSIONS
        """
    } else {
        // URL download — existing logic
        def url_basename = source.split('\\?')[0].split('/')[-1]
        def ext = url_basename.endsWith('.fna.gz') ? '.fna.gz' :
                  url_basename.endsWith('.fa.gz')  ? '.fa.gz' :
                  url_basename.endsWith('.fasta.gz') ? '.fasta.gz' :
                  url_basename.endsWith('.fna') ? '.fna' :
                  url_basename.endsWith('.fasta') ? '.fasta' :
                  url_basename.endsWith('.fa') ? '.fa' :
                  '.fasta'  // default for API URLs with no extension
        """
        wget -q -O "${meta.id}_downloaded${ext}" "${source}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            wget: \$(wget --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+' || echo "1.21")
        END_VERSIONS
        """
    }
}
