process COGEQC {
    tag "Orthogroup Summary"
    label 'process_medium'

    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/cogeqc-1.2.0_r-4.2.2':
        '' }"

    publishDir(
        path: "${params.outdir}/orthogroup-summaries",
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path orthofinder_outdir
    path prot_annotations // Base filepath to where protein annotations are stored

    output:
    path "*-cogeqc-summary.tsv" , emit: og_summary
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    """
    # Assess orthogroups inferred using each inflation parameter, summarizing
    # how well they group proteins with the same domains together, as well as
    # other summary stats like number of ogs with >= 4 species, per-species
    # gene count per-og, etc.
    Rscript $projectDir/bin/cogeqc-summarize-ogs.R ${orthofinder_outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cogeqc: \$( cat version.txt | head -n1 | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}
