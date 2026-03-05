// Rename FASTA files so their filename matches the normalized species name.
// OrthoFinder uses filenames as species identifiers, so this ensures
// all downstream tip labels use the hyphenated species name convention
// (e.g. Homo-sapiens_ProteinID).
process RENAME_FASTAS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: renamed

    script:
    """
    cp ${fasta} ${meta.id}.fa
    """
}
