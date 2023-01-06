process ANNOTATE_UNIPROT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/r-4.2.2_uniprot.ws:2.38.0':
        '' }"

    publishDir(
        path: "${params.outdir}/protein-annotations",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(annots_to_download), val(meta), path(fasta)

    output:
    path "*accessions.txt"  , emit: accessions
    path "*" , emit: all_annotations
    path "*cogeqc-annotations.tsv" , emit: cogeqc_annotations
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    def spp = "${meta.id}"
    def is_uniprot = "${meta.uniprot}"
    """
    # Only annotate species for which protein IDs are found in UniProt (i.e.
    # proteomes come from UniProt).
    # Check below - if from uniprot, go ahead and annotate, otherwise skip the species.
    if [ "$is_uniprot" == "true" ]; then
        # Pull out the sequence names, strip trailing info, and remove spp name.
        grep ">" $fasta | cut -d" " -f1 | cut -d":" -f2 > $spp-protein-accessions.txt

        # Now run the script to pull down annotations for the protein accessions in this species.
        # This R script uses the UniProt.ws bioconducter package to accomplish this.
        # NOTE: The script is packaged in the bin/ subdirectory of this workflow.

        UniProt-Protein-Annotation-NF.R $annots_to_download $spp ${spp}-protein-accessions.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | sed "s/R version //g" | sed "s/ (.*//g" )
        UniProt.ws: \$( cat version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}
