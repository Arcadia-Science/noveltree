process GENERAX_PER_SPECIES {
    tag "GeneRax: Per-species rates"
    label 'process_generax_per_species'
    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems. 
    container "${ workflow.containerEngine == 'docker' ?
        'arcadiascience/generax_19604b7:0.0.1': '' }"

    publishDir(
        path: "${params.outdir}/generax/per_species_rates",
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
        --families generax_orthogroup.families \\
        --per-species-rates \\
        --prefix GeneRax \\
        $args

    # And move the results into the current working directory
    mv GeneRax/* .
    rm -r GeneRax
    rm -r results
    """
}
