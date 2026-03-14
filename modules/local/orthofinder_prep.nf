process ORTHOFINDER_PREP {
    tag "Prepping data for OrthoFinder"
    label 'process_low'

    container 'arcadiascience/orthofinder_2.5.4:1.0.0'

    storeDir "${params.outdir}/orthofinder_prep"

    input:
    file(fasta)
    val output_directory

    output:
    path "**.dmnd"           , emit: diamonds
    path "**.fa"             , emit: fastas
    path "**SequenceIDs.txt" , emit: seqIDs
    path "**SpeciesIDs.txt"  , emit: sppIDs

    script:
    """
    # The fasta directory depends on whether we're running the mcl testing or not.
    orthofinder \\
        -f ./ \\
        -t ${task.cpus} \\
        -op > tmp

    mkdir ${output_directory} && mv OrthoFinder/ ${output_directory}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder --versions | head -n2 | tail -n1 | sed "s/OrthoFinder version //g" | sed "s/ Copyright (C) 2014 David Emms//g")
    END_VERSIONS
    """
}
