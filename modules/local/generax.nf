process GENERAX {
    tag "GeneRax"
    label 'process_highthread'

    conda (params.enable_conda ? "bioconda::generax==2.0.4--h19e7193_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/generax:2.0.4--h19e7193_0':
        'quay.io/biocontainers/generax:2.0.4--h19e7193_0' }"
        
    publishDir(
        path: "${params.outdir}/generax",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path init_species_tree // Filepath to the SpeciesRax species tree
    path generax_map       // Filepath to the generax gene-species map file
    path gene_trees        // Filepaths to the starting gene trees
    path alignments        // Filepaths to the gene family alignments
    path families          // Filepath to the families file

    output:
    path "*" , emit: results

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''

    """
    mpiexec -np ${task.cpus} --allow-run-as-root generax \
    --species-tree $init_species_tree \
    --families $families \
    --rec-model UndatedDTL \
    --prune-species-tree \
    --per-family-rates \
    --strategy SPR \
    --prefix GeneRax
    """
}

