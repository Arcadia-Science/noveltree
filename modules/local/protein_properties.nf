process PROTEIN_PROPERTIES {
    tag "${meta.og}"

    cpus 1
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def estimated_gb = Math.max(4L, (long)(n * L * 32L / (1024L * 1024L * 1024L)) + 1L)
        def capped_gb = (int) Math.min(estimated_gb, 16L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/protein_properties:1.0.0'

    storeDir "${params.outdir}/protein_properties"

    input:
    tuple val(meta), path(original_fasta), path(cleaned_msa)

    output:
    tuple val(meta), path("aa-summary-stats/per-family-summaries/aa-physical-properties/${meta.og}_summary_statistics.csv"), emit: summary_stats
    tuple val(meta), path("aa-summary-stats/per-family-summaries/aa-physical-property-sds/${meta.og}_summary_statistics_sd.csv"), emit: summary_stats_sd
    tuple val(meta), path("aa-summary-stats/per-family-summaries/aa-physical-property-autocorr/${meta.og}_summary_statistics_autocorr.csv"), emit: summary_stats_autocorr
    path "aa-summary-stats/per-family-summaries/aa-counts/${meta.og}_aa_composition_counts.csv"           , emit: aa_counts
    path "aa-summary-stats/per-family-summaries/aa-proportions/${meta.og}_aa_composition_percentages.csv" , emit: aa_proportions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract survivor protein IDs from the cleaned MSA headers
    grep ">" ${cleaned_msa} | sed 's/>//' | cut -d' ' -f1 > survivors.txt

    # Filter original FASTA to only proteins that survived alignment trimming
    awk 'BEGIN{while((getline line < "survivors.txt") > 0) ids[line]=1}
         /^>/{p=ids[substr(\$1,2)]} p' ${original_fasta} > ${meta.og}_filtered.fa

    # Run property calculation on filtered full-length sequences
    genefam_aa_summaries.py ${meta.og}_filtered.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
