process PHYSICOCHEMICAL_PROPS {
    tag "${meta.og}"
    label "process_high"

    container 'arcadiascience/physicochemical_props:1.0.0'

    input:
    tuple val(meta), path(msa_file)

    output:
    tuple val(meta), path("aa-summary-stats/per-family-summaries/aa-physical-properties/${meta.og}_summary_statistics.csv"), emit: summary_stats
    path "aa-summary-stats/per-family-summaries/aa-counts/${meta.og}_aa_composition_counts.csv"           , emit: aa_counts
    path "aa-summary-stats/per-family-summaries/aa-proportions/${meta.og}_aa_composition_percentages.csv" , emit: aa_proportions
    path "versions.yml"                                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Run the physicochemical properties calculation script
    # Script processes a single MSA file and outputs with gene family name in filename
    genefam_aa_summaries.py ${msa_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
