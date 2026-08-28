// Normalize FASTA files and sequence identifiers once at the input boundary.
// OrthoFinder uses filenames as species identifiers, so filenames are set
// to the normalized species name. The emitted mapping is the authoritative
// source of gene-to-species identity for downstream reconciliation.
process RENAME_FASTAS {
    tag "${meta.id}"
    label 'process_single'
    container 'arcadiascience/preprocess_proteomes:1.1.0'
    maxRetries 0
    errorStrategy 'terminate'

    storeDir "${params.outdir}/preprocessing/final_proteomes"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: renamed
    tuple val(meta), path("${meta.id}.protein_map.tsv"), emit: protein_map
    tuple val(meta), path("${meta.id}.normalization_qc.json"), emit: normalization_qc

    script:
    """
    normalize_proteome_fasta.py \
        --input ${fasta} \
        --species '${meta.id}' \
        --output ${meta.id}.fa \
        --mapping ${meta.id}.protein_map.tsv \
        --qc ${meta.id}.normalization_qc.json
    """
}
