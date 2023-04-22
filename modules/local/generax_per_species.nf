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

    input:
    file species_tree      // Filepath to the SpeciesRax species tree
    file generax_map       // Filepath to the generax gene-species map file
    file gene_trees        // Filepaths to the starting gene trees
    file alignments        // Filepaths to the gene family alignments
    file families          // Filepath to the families file

    output:
    path "*"                        , emit: results

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

    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $species_tree \\
        --families $families \\
        --per-species-rates \\
        --prefix GeneRax \\
        $args

    # And move the results into the current working directory
    mv GeneRax/* .
    rm -r GeneRax
    rm -r results
    """
}
