process ORTHOFINDER_PREP {
    tag "Prepping data for OrthoFinder"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.4--hdfd78af_0' :
        'quay.io/biocontainers/orthofinder:2.5.4--hdfd78af_0' }"

    //publishDir(
    //    path: "${params.outdir}/orthofinder/",
    //    mode: 'copy',
    //    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    //)
    
    input:
    tuple val(meta), path(fasta)    // samplesheet with filepaths
    val mcl_test
    
    output:
    path "*.dmnd", emit: dmd
    path "*.fa", emit: fa
    path "SequenceIDs.txt", emit: seqIDs
    path "SpeciesIDs.txt", emit: sppIDs
    path "versions.yml" , emit: versions

    script:
    def fas = path(fasta)
    def testing_mcl = mcl_test.equals('true') ? "${mcl_test}" : "false"
    """
    # Prepare the bespoke orthofinder directory structure and stripped down files
    # spit the output commands to a tmp file - we don't care about this
    input_dir=\$(cat ./)
    orthofinder \\
        -f \$input_dir \\
        -t ${task.cpus} \\
        -op > tmp
        
    # Create the publishdir if it's not made yet
    mkdir -p ../../../${params.outdir}/orthofinder/
    
    # Copy all the other orthofinder scraps (species and sequence IDs, etc) to
    # here - needed for downstream interfacing with orthofinder. 
    # Specific outdir will depend on whether we're testing mcl or not. 
    if [ "$testing_mcl" == "true" ]; then
        outDir="mcl_testing"
    else
        outDir="full_analysis"
    fi
    
    mkdir -p ../../../${params.outdir}/orthofinder/\$outDir/
    mkdir -p ../../../${params.outdir}/orthofinder/\$outDir/data
    dataDir="../../../${params.outdir}/orthofinder/\$outDir/data"
    cp \$input_dir/OrthoFinder/Results*/WorkingDirectory/* \$dataDir/
    cp \$input_dir/OrthoFinder/Results*/WorkingDirectory/* .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder --versions | head -n2 | tail -n1 | sed "s/OrthoFinder version //g" | sed "s/ Copyright (C) 2014 David Emms//g")
    END_VERSIONS
    """
}
