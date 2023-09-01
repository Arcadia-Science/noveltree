process BUSCO {
    // Modified from nf-core to:
    // 1) specify required "mode" parameter
    // 2) allow scale-dependent (from meta) specification of lineage dataset
    tag "$meta.id"
    label 'process_low_cpu'

    conda (params.enable_conda ? "bioconda::busco=5.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.4.3--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.4.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path('tmp_input/*')
    val lineage_scale                     // Required: Taxonomically shallow- or broad-scale busco database to compare against? Determines busco lineage to use from samplesheet.
    path busco_lineages_path              // Recommended: path to busco lineages - downloads if not set
    path config_file                      // Optional:    busco configuration file

    output:
    tuple val(meta), path("*_busco.batch_summary.txt") , emit: batch_summary, optional: true
    tuple val(meta), path("short_summary.*.txt")       , emit: short_summaries_txt, optional: true
    tuple val(meta), path("short_summary.*.json")      , emit: short_summaries_json, optional: true
    tuple val(meta), path("*_busco.tar.gz")            , emit: busco_dir, optional: true
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args              = task.ext.args ?: ''
    def prefix            = lineage_scale.equals('shallow') ? "${meta.id}_${meta.shallow_db}" : "${meta.id}_${meta.broad_db}"
    def busco_config      = config_file ? "--config $config_file" : ''
    def busco_lineage     = lineage_scale.equals('shallow') ? "--lineage_dataset ${meta.shallow_db}" : "--lineage_dataset ${meta.broad_db}"
    def busco_lineage_dir = busco_lineages_path ? "--offline --download_path ${busco_lineages_path}" : ''
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    busco \\
        --cpu ${task.cpus} \\
        --in "\$INPUT_SEQS" \\
        --out ${prefix}_busco \\
        --mode ${meta.mode} \\
        $busco_lineage \\
        $busco_lineage_dir \\
        $busco_config \\
        $args

    # clean up
    rm -rf "\$INPUT_SEQS"

    # Move files to avoid staging/publishing issues
    mv ${prefix}_busco/batch_summary.txt ${prefix}_busco.batch_summary.txt
    mv ${prefix}_busco/*/short_summary.*.{json,txt} . || echo "Short summaries were not available: No genes were found."

    # Compress and then remove the busco dir, as these contain many largely unnecessary files:
    tar -czvf ${prefix}_busco.tar.gz ${prefix}_busco/
    rm -r ${prefix}_busco/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
