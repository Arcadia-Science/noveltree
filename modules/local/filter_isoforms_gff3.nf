process FILTER_ISOFORMS_GFF3 {
    tag "${meta.id}"
    label 'process_single'
    container 'arcadiascience/preprocess_proteomes:1.1.0'
    errorStrategy 'terminate'
    maxRetries 0

    storeDir "${params.outdir}/preprocessing/isofiltered"

    input:
    tuple val(meta), path(fasta), path(annotation, stageAs: 'annotation.gff3')

    output:
    tuple val(meta), path("${meta.id}_isofiltered.fasta"), emit: filtered
    tuple val(meta), path("${meta.id}_isoform_selection.tsv"), emit: ledger
    tuple val(meta), path("${meta.id}_isoform_qc.json"), emit: qc

    script:
    """
    filter_isoforms_gff3.py --fasta "${fasta}" --gff3 "${annotation}" \\
        --species '${meta.id}' --output '${meta.id}_isofiltered.fasta' \\
        --ledger '${meta.id}_isoform_selection.tsv' --qc '${meta.id}_isoform_qc.json'
    """
}
