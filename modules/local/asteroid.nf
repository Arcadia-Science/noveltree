process ASTEROID {
    tag "Asteroid"
    label 'process_highcpu'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/asteroid_3aae117:0.0.1':
        '' }"

    publishDir(
        path: "${params.outdir}/asteroid",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path treefile     // Filepath to the asteroid treefile (all newick gene trees)
    path asteroid_map // Filepath to the asteroid gene-species map

    output:
    path "*bestTree.newick" , emit: spp_tree
    path "*allTrees.newick" , emit: all_asteroid_trees
    path "*scores.txt"      , emit: asteroid_scores
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    """
    # Run asteroid using the multithreaded mpi version
    mpiexec -np ${task.cpus} --allow-run-as-root --use-hwthread-cpus \
    asteroid \
    -i $treefile \
    -m $asteroid_map \
    -p asteroid \
    $args

    # Version is hardcoded for now (asteroid doesn't output this currently)
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        asteroid: 1.0
    END_VERSIONS
    """
}
