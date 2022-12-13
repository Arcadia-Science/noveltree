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
    path fasta_dir // directory containing fasta files
    val out_dir   // name of directory where prepped files will be stored
    val fasta_list_fname
    val dmnd_list_fname
    
    output:
    path fasta_list_fname, emit: fa
    path dmnd_list_fname, emit: dmd
    path "*.dmnd", emit: dmnd
    path "*.fa", emit: fastas
    path "SequenceIDs.txt", emit: seqIDs
    path "SpeciesIDs.txt", emit: sppIDs
    path "versions.yml" , emit: versions

    script:
    """
    # The fasta directroy depends on whether we're running the mcl testing or not.  
    fastaDir=\$(cat $fasta_dir)
    orthofinder \\
        -f \$fastaDir \\
        -t ${task.cpus} \\
        -op > tmp
        
    # Move the output to a more sensible directory
    mkdir -p ${params.outdir}/orthofinder/$out_dir
    cp \${fastaDir}/OrthoFinder/Results*/WorkingDirectory/* ${params.outdir}/orthofinder/$out_dir
    cp \${fastaDir}/OrthoFinder/Results*/WorkingDirectory/* .
    rm -r \${fastaDir}/OrthoFinder/
    
    # output the paths to fasta files and dmnd dbs
    ls ${params.outdir}/orthofinder/$out_dir/*fa > $fasta_list_fname
    ls ${params.outdir}/orthofinder/$out_dir/*dmnd > $dmnd_list_fname

    #mkdir -p ../../../${params.outdir}/orthofinder/\$outDir/
    #mkdir -p ../../../${params.outdir}/orthofinder/\$outDir/data
    #dataDir="../../../${params.outdir}/orthofinder/\$outDir/data"
    #cp \$input_dir/OrthoFinder/Results*/WorkingDirectory/* \$dataDir/
    #cp \$input_dir/OrthoFinder/Results*/WorkingDirectory/* .
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder --versions | head -n2 | tail -n1 | sed "s/OrthoFinder version //g" | sed "s/ Copyright (C) 2014 David Emms//g")
    END_VERSIONS
    """
}
