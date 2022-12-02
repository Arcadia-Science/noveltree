process ANNOTATE_UNIPROT {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/r-4.2.2_uniprot.ws:2.38.0':
        '' }"
        
    publishDir(
        path: "${params.outdir}/protein-annotations",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)    // path('tmp_input/*')

    output:
    path "accessions/*"  , emit: accessions
    path "annotations/*" , emit: annotations
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
        # First pull out the protein IDs and write these out to file. 
        # Make a directory to store these.
        mkdir accessions
        
        # Pull out the sequence names, strip trailing info, and remove spp name. 
        grep ">" $fasta | cut -d" " -f1 | cut -d":" -f2 > accessions/$spp-protein-accessions.txt
        
        # Now run the script to pull down annotations for the protein accessions in this species. 
        # This R script uses the UniProt.ws bioconducter package to accomplish this. 
        # NOTE: The script is packaged in the bin/ subdirectory of this workflow. 
        
        Rscript $projectDir/bin/UniProt-Protein-Annotation-NF.R $spp accessions/${spp}-protein-accessions.txt
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$( R --version | head -n1 | sed "s/R version //g" | sed "s/ (.*//g" )
        UniProt.ws: \$( cat version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}
