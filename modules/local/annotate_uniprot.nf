process ANNOTATE_UNIPROT {
    tag "$meta.id"
    label 'process_medium'

    container 'arcadiascience/bioservices_1.10.0:1.0.0'

    publishDir(
        path: "${params.outdir}/protein_annotations",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), file(fasta)

    output:
    path "*accessions.txt"         , emit: accessions
    path "*"                       , emit: all_annotations
    path "*cogeqc_annotations.tsv" , emit: cogeqc_annotations
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def spp         = "${meta.id}"
    def is_uniprot  = "${meta.uniprot}"
    def project_dir = "${projectDir}"
    """
    # Only annotate species for which protein IDs are found in UniProt (i.e.
    # proteomes come from UniProt).
    # Check below - if from uniprot, go ahead and annotate, otherwise skip the species.
    if [ "$is_uniprot" == "true" ]; then
        # Pull out the sequence names, strip trailing info, and remove spp name.
        # Handle both colon-delimited (>Species:Accession) and pipe-delimited (>Species|Accession|Entry) formats
        grep ">" $fasta | cut -d" " -f1 | awk -F'[:|]' '{print \$2}' > ${spp}_protein_accessions.txt

        # Retrieve InterPro annotations from UniProt REST API.
        protein_annotation.py $spp ${spp}_protein_accessions.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$( python --version | sed "s/Python //g" | sed "s/ (.*//g" )
        requests: \$(python -c "import requests; print(requests.__version__)")
    END_VERSIONS
    """
}
