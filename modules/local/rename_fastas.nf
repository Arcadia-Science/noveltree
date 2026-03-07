// Rename FASTA files and normalize sequence headers.
// OrthoFinder uses filenames as species identifiers, so filenames are set
// to the normalized species name. Headers are rewritten as
// >{species_name}_{sanitized_id} where underscores in the protein ID are
// replaced with hyphens. This ensures all downstream tools can extract the
// species name by stripping everything after the last underscore.
process RENAME_FASTAS {
    tag "${meta.id}"
    label 'process_single'
    container 'arcadiascience/rbase_4.2.2:1.0.0'
    maxRetries 0
    errorStrategy 'terminate'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: renamed

    script:
    """
    awk -v species="${meta.id}" '
    /^>/ {
        # Extract first word (sequence ID)
        id = \$1
        sub(/^>/, "", id)

        # If already prefixed with species name, strip it
        prefix = species "_"
        if (index(id, prefix) == 1) {
            prot = substr(id, length(prefix) + 1)
        } else {
            prot = id
        }

        # Replace underscores with hyphens in protein part
        gsub(/_/, "-", prot)

        print ">" species "_" prot
        next
    }
    { print }
    ' ${fasta} > ${meta.id}.fa
    """
}
