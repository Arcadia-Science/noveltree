process PROCESS_ZOOGLE_TABLE {
    tag "Process zoogle table"
    label "process_medium"

    container "${
        (workflow.containerEngine == 'docker') || (workflow.containerEngine == 'singularity') ?
        'quay.io/biocontainers/python:3.9--1':''
    }"

    publishDir(
        path: "${params.outdir}/zoogle_processed",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    )

    input:
    path aggregated_table  // all_protein_comparisons.tsv.gz from AGGREGATE_PHYLO_DIST
    path hgnc_processed    // Processed HGNC dataset from DOWNLOAD_HGNC

    output:
    path "zoogle_final_table.tsv.gz"  , emit: final_table
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Install required Python packages
    pip install --no-cache-dir pandas click tqdm

    # Add zoogle-specific columns to the aggregated table
    # This uses the exact same processing functions from the 2025-zoogle repository
    python ${projectDir}/bin/add_zoogle_columns.py \\
        --input ${aggregated_table} \\
        --hgnc ${hgnc_processed} \\
        --output zoogle_final_table.tsv.gz

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        click: \$(python -c "import click; print(click.__version__)")
        tqdm: \$(python -c "import tqdm; print(tqdm.__version__)")
    END_VERSIONS
    """
}