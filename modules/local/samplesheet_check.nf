process SAMPLESHEET_CHECK {
    tag "$complete_samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"
    
    // TODO: FOR SOME REASON THIS MODULE OUTPUTS TO "PIPELINE_INFO" BUT THE MCL_SAMPLESHEET_CHECK OUTPUTS TO ITS OWN DIRECTORY. 
    // SPECIFYNG PUBLISHDIR BELOW DOESN'T SEEM TO FIX IT?
    //publishDir(
    //    path: "${params.outdir}/complete_dataset",
    //    mode: 'copy',
    //    saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    //)

    input:
    path complete_samplesheet // Samplesheet formatted as described in the README
    path data_dir
    
    output:
    path complete_samplesheet    , emit: complete_csv     // Original samplesheet
    path 'mcl_test_samplesheet.csv' , emit: mcl_test_csv   // Per-species proteomes
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This python script is bundled with the pipeline, in phylorthology/bin/

    """
    check_samplesheet.py \\
        $complete_samplesheet \\
        complete_samplesheet.valid.csv
    
    mkdir -p ${params.outdir}/complete_dataset
    mkdir -p ${params.outdir}/mcl_test_dataset
    mv $data_dir/*.fasta ${params.outdir}/complete_dataset
    
    # Create the mcl testing datasheet and copy over to its directory
    head -n1 $complete_samplesheet > mcl_test_samplesheet.csv
    awk -F, '\$8 == "true"' $complete_samplesheet >> mcl_test_samplesheet.csv
    
    for spp in \$(tail -n+2 mcl_test_samplesheet.csv | cut -f1 -d",")
    do
        cp ${params.outdir}/complete_dataset/\${spp}* ${params.outdir}/mcl_test_dataset/
    done 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
