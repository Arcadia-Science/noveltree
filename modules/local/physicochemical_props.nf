process PHYSICOCHEMICAL_PROPS {
    tag "Physicochemical Properties"
    label "process_high"

    container 'arcadiascience/physicochemical_props:1.0.0'

    input:
    path msa_files

    output:
    path "aa-summary-stats/across-family-summaries/*.csv"   , emit: across_family_summaries
    path "aa-summary-stats/per-family-summaries/**/*.csv"   , emit: per_family_summaries
    path "aa-summary-stats/"                                , emit: all_outputs
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Run the physicochemical properties calculation script
    # Nextflow stages all input files into the work directory
    # The script expects a directory path containing the MSA files
    genefam_aa_summaries.py . --threads ${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
