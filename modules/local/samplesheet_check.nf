process SAMPLESHEET_CHECK {
    tag "$complete_samplesheet"
    errorStrategy 'terminate'
    maxRetries 0

    container 'arcadiascience/python_3.9'

    input:
    path complete_samplesheet // Samplesheet formatted as described in the README

    output:
    path '*.csv'        , emit: csv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This python script is bundled with the pipeline, in bin/
    """
    check_samplesheet.py \\
        $complete_samplesheet \\
        complete_samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
