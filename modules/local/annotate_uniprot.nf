process ANNOTATE_UNIPROT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/bioservices_1.10.0:1.0.0':
        '' }"

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
    def args               = task.ext.args ?: ''
    def spp                = "${meta.id}"
    def is_uniprot         = "${meta.uniprot}"
    def project_dir        = "${projectDir}"
    def sppid_protid_delim = "${params.sppid_protid_delim}"
    """
    # Only annotate species for which protein IDs are found in UniProt (i.e.
    # proteomes come from UniProt).
    # Check below - if from uniprot, go ahead and annotate, otherwise skip the species.
    if [ "$is_uniprot" == "true" ]; then
        # Pull out the sequence names, strip trailing info, and remove spp name.
        grep ">" $fasta | cut -d" " -f1 | cut -d"$sppid_protid_delim" -f2 > ${spp}_protein_accessions.txt
        
        # If the data was pre-processed using our snakemake workflow, the uniprot
        # accessions will actually be nested between to pipes ("|") - for example:
        # "tr|UNIPROTID|UNIPROTID-SppAbbrev"
        # Make sure that the file containing protein accessions have pulled these out
        cut -f2 -d"|" ${spp}_protein_accessions.txt > tmp && mv tmp ${spp}_protein_accessions.txt
        
        # Now run the script to pull down annotations for the protein accessions in this species.
        # This Python script uses the bioservices python package to accomplish this.
        # NOTE: The script is packaged in the bin/ subdirectory of this workflow.
        protein_annotation.py $spp ${spp}_protein_accessions.txt

        # Organize results so that the cogeqc annotations are in the current
        # directory, and all others are moved into a single directory for the
        # species
        mkdir $spp
        for f in \$(ls *.tsv | grep -v "cogeqc")
        do
            mv \$f ${spp}/
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$( python --version | sed "s/Python //g" | sed "s/ (.*//g" )
        bioservices: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('bioservices').version)")
    END_VERSIONS
    """
}
