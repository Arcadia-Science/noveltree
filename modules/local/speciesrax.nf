process SPECIESRAX {
    tag "SpeciesRax"
    //label 'process_highthread' // Possible specification for full analysis
    label 'process_medium' // Used for debugging

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
    path init_species_tree // Filepath to the starting species tree
    path generax_map       // Filepath to the generax gene-species map file
    path gene_trees        // Filepaths to the starting gene trees
    path alignments        // Filepaths to the gene family alignments
    path families          // Filepath to the families file

    output:
    path "*" ,                                    emit: results
    path "speciesrax_final_species_tree.newick" , emit: speciesrax_tree
    path "versions.yml" ,                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''

    """
    mpiexec -np ${task.cpus} --allow-run-as-root generax \
    --species-tree $init_species_tree \
    --families $families \
    --per-family-rates \
    --rec-model UndatedDTL \
    --prune-species-tree \
    --si-strategy HYBRID \
    --si-quartet-support \
    --si-estimate-bl \
    --strategy SPR \
    --prefix SpeciesRax
    
    # Copy the final starting tree to the current working directory - 
    # this will be used as the starting tree for generax application to the 
    # remaining gene families. 
    cp SpeciesRax/species_trees/inferred_species_tree.newick ./speciesrax_final_species_tree.newick

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$( generax | head -n1 | sed "s/.*GeneRax //g" )
    END_VERSIONS
    """
}

