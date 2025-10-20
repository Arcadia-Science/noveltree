process DOWNLOAD_HGNC {
    tag "Download and process HGNC"
    label "process_low"

    container "${
        (workflow.containerEngine == 'docker') || (workflow.containerEngine == 'singularity') ?
        'quay.io/biocontainers/python:3.9--1':''
    }"

    publishDir(
        path: "${params.outdir}/hgnc",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    )

    output:
    path "hgnc_processed.tsv"  , emit: hgnc_processed
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install required Python packages
    pip install --no-cache-dir pandas click

    # Download the HGNC dataset
    python ${projectDir}/bin/hgnc_dataset.py download \\
        --output-filepath hgnc_raw.tsv

    # Process the HGNC dataset
    python ${projectDir}/bin/hgnc_dataset.py process \\
        --input-filepath hgnc_raw.tsv \\
        --output-filepath hgnc_processed.tsv

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        click: \$(python -c "import click; print(click.__version__)")
    END_VERSIONS
    """
}