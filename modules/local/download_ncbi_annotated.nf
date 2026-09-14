// Fetch both files from the SAME assembly package. A protein-only cached
// download cannot satisfy this process's required annotation/metadata outputs.
process DOWNLOAD_NCBI_ANNOTATED {
    tag "${meta.id}"
    label 'process_single'
    container 'arcadiascience/preprocess_proteomes:1.1.0'
    storeDir "${params.outdir}/downloads"

    input:
    tuple val(meta), val(source)

    output:
    tuple val(meta), path("${meta.id}_downloaded.fasta"), emit: downloaded
    tuple val(meta.id), path("${meta.id}_annotation.gff3"), emit: annotation
    tuple val(meta.id), path("${meta.id}_download_metadata.json"), emit: metadata

    script:
    """
    datasets download genome accession '${source}' --include protein,gff3 --filename data.zip
    unzip -q data.zip -d ncbi_data
    extract_ncbi_proteome.py --directory ncbi_data --accession '${source}' \\
        --protein-output '${meta.id}_downloaded.fasta' \\
        --annotation-output '${meta.id}_annotation.gff3' \\
        --metadata-output '${meta.id}_download_metadata.json'
    """
}
