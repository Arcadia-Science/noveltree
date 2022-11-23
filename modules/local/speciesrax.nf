process SPECIESRAX {
    tag "SpeciesRax"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::generax==2.0.4--h19e7193_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/generax:2.0.4--h19e7193_0':
        'quay.io/biocontainers/generax:2.0.4--h19e7193_0' }"
        
    publishDir(
        path: "${params.outdir}/speciesrax",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path init_species_tree   // Filepath to the starting species tree
    path generax_map // Filepath to the generax gene-species map file
    path gene_trees // Filepaths to the starting gene trees
    path alignments // Filepaths to the gene family alignments
    path families // Filepath to the families file

    output:
    path "*" , emit: results
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    mpiexec -np 32 --allow-run-as-root generax --families $families --strategy SPR \
    --si-strategy HYBRID --species-tree $init_species_tree --rec-model UndatedDTL \
    --per-family-rates --prune-species-tree --si-estimate-bl \
    --si-spr-radius 5 --max-spr-radius 5 --si-quartet-support \
    --prefix SpeciesRax

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$( generax | head -n1 | sed "s/.*GeneRax //g" )
    END_VERSIONS
    """
}

