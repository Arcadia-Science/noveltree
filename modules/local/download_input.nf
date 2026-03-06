process DOWNLOAD_INPUT {
    tag "${meta.id}"
    label 'process_single'

    container 'arcadiascience/preprocess_proteomes:1.0.0'

    input:
    tuple val(meta), val(source_url)

    output:
    tuple val(meta), path("${meta.id}_downloaded.*"), emit: downloaded
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Determine output extension from URL
    def url_basename = source_url.split('\\?')[0].split('/')[-1]
    def ext = url_basename.endsWith('.fna.gz') ? '.fna.gz' :
              url_basename.endsWith('.fa.gz')  ? '.fa.gz' :
              url_basename.endsWith('.fasta.gz') ? '.fasta.gz' :
              url_basename.endsWith('.fna') ? '.fna' :
              url_basename.endsWith('.fasta') ? '.fasta' :
              url_basename.endsWith('.fa') ? '.fa' :
              '.fasta'  // default for API URLs with no extension
    """
    wget -q -O "${meta.id}_downloaded${ext}" "${source_url}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version 2>&1 | head -1 | grep -oP '\\d+\\.\\d+' || echo "1.21")
    END_VERSIONS
    """
}
