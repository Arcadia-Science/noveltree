process ORTHOFINDER_PREP {
    tag "Prepping data for OrthoFinder"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.3--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.3--hdfd78af_0' }"

    publishDir(
        path: "${params.outdir}/orthofinder/WorkingDirectory",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
    input:
    val fasta_dir
    
    output:
    path "*.dmnd", emit: dmd
    path "*.fa", emit: fa
    path "SequenceIDs.txt", emit: seqIDs
    path "SpeciesIDs.txt", emit: sppIDs
    path "versions.yml" , emit: versions

    script:
    """
    # Prepare the bespoke orthofinder directory structure and stripped down files
    # spit the output commands to a tmp file - we don't care about this
    orthofinder \\
        -f $fasta_dir \\
        -t $task.cpus \\
        -op > tmp
        
    # Copy all the other orthofinder scraps (species and sequence IDs, etc) to
    # here - needed for downstream interfacing with orthofinder. 
    cp ${fasta_dir}/OrthoFinder/Results*/WorkingDirectory/* .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder --versions | head -n2 | tail -n1 | sed "s/OrthoFinder version //g" | sed "s/ Copyright (C) 2014 David Emms//g")
    END_VERSIONS
    """
}
