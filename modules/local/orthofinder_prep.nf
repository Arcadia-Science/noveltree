process ORTHOFINDER_PREP {
    tag "Prepping data for OrthoFinder"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.4--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.4--hdfd78af_0' }"

    input:
    path(fasta), stageAs: "fasta/"

    output:
    // path fasta_list_fname, emit: fa
    // path dmnd_list_fname, emit: dmd
    path "**.dmnd", emit: dmnd
    path "**.fa", emit: fastas
    // path "SequenceIDs.txt", emit: seqIDs
    // path "SpeciesIDs.txt", emit: sppIDs
    path "versions.yml" , emit: versions

    script:
    """
    # The fasta directroy depends on whether we're running the mcl testing or not.
    orthofinder \\
        -f fasta/ \\
        -t ${task.cpus} \\
        -op > tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder --versions | head -n2 | tail -n1 | sed "s/OrthoFinder version //g" | sed "s/ Copyright (C) 2014 David Emms//g")
    END_VERSIONS
    """
}
