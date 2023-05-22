process GENERAX_PER_FAMILY {
    tag "$meta.og"
    label 'process_generax'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems. 
    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/generax_19604b7:0.0.1': '' }"

    publishDir(
        path: "${params.outdir}/generax/per_family_rates",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input: // Input is a single large tuple with paths to map-links, tree files, alignments, and the species tree
    tuple val(meta), file(map_link), file(gene_tree), file(alignment), file(species_tree)

    output:
    path "*"                                         , emit: results
    tuple val(meta), path("**_reconciled_gft.newick"), emit: generax_per_fam_gfts

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    def og      = "${meta.og}"
    """
    # Recode selenocysteine as a gap character:
    # RAxML-NG (used under the hood by SpeciesRax and
    # GeneRax) cannot handle these. Even if rare,
    # their inclusion leads a number of gene families
    # to be excluded from analyses.
    sed -E -i '/>/!s/U/-/g' *.fa

    # Do the same for Pyrrolysine
    sed -E -i '/>/!s/O/-/g' *.fa

    # Populate the family file for this gene family for the 
    # analysis with GeneRax
    # We will be using LG+G4+F for all gene families
    echo "[FAMILIES]" > ${og}.family
    echo "- ${og}" >> ${og}.family
    echo "starting_gene_tree = ${gene_tree}" >> ${og}.family
    echo "mapping = ${og}_map.link" >> ${og}.family
    echo "alignment = $alignment" >> ${og}.family
    echo "subst_model = LG+G4+F" >> ${og}.family

    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $species_tree \\
        --families ${og}.family \\
        --per-family-rates \\
        --prefix $og \\
        --prune-species-tree \\
        --reconciliation-samples 100 \\
        $args

    # Clean up
    rm -r $og/gene_optimization_*

    # Rename the inferred reconciled gene trees to be named after their corresponding orthogroup
    mv $og/results/$og/geneTree.newick $og/results/$og/${og}_reconciled_gft.newick
    
    # And move the reconciliation transfer samples into a subdirectory, archive, and compress.
    mkdir $og/reconciliations/reconciliation_transfer_samples/
    mv $og/reconciliations/*_*_transfers.txt $og/reconciliations/reconciliation_transfer_samples/
    tar -czvf $og/reconciliations/reconciliation_transfer_samples.tar.gz $og/reconciliations/reconciliation_transfer_samples/
    """
}
