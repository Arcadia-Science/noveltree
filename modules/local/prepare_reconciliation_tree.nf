process PREPARE_RECONCILIATION_TREE {
    tag "${meta.og}"
    label 'process_single'

    container 'arcadiascience/python_3.9'
    storeDir "${params.outdir}/gene_family_trees/reconciliation_ready"

    input:
    tuple val(meta), path(gene_tree), path(mapping), path(alignment)

    output:
    tuple val(meta), path("${meta.og}_reconciliation_ready.newick"), emit: tree
    tuple val(meta), path("${meta.og}_reconciliation_qc.tsv"), emit: qc

    script:
    """
    prepare_reconciliation_tree.py \
        --orthogroup '${meta.og}' \
        --tree ${gene_tree} \
        --alignment ${alignment} \
        --mapping ${mapping} \
        --output ${meta.og}_reconciliation_ready.newick \
        --qc ${meta.og}_reconciliation_qc.tsv
    """
}
