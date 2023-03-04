process CIALIGN {
    tag "$fasta"
    label 'process_lowcpu'

    // container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/cialign_1.1.0:0.0.1' :
    //     '' }"
    container "${ workflow.containerEngine == 'docker' ? 'austinhpatton123/cialign_1.1.0:0.0.1' :
        '' }"

    publishDir(
        path: "${params.outdir}/cialign_cleaned_msas",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(fasta)   // Filepaths to the MSAs

    output:
    path("*_cialign.fa") , emit: trimmed_msas
    path "*"             , emit: results
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Get the name of the orthogroup we are processing
    prefix=\$(echo $fasta | cut -f1 -d "_")

    # Clean up the MSAs for each orthogroup containing at least 4 species.
    CIAlign \
        --infile ${fasta} \
        --outfile_stem="\${prefix}" \
        $args
    
    # Rename output so it is clear we trimmed with CIAlign
    mv \${prefix}_cleaned.fasta \${prefix}_cialign.fa
    
    # And move the "removed.txt" files indicating which sites were removed 
    # from each MSA to a separate directory
    mkdir removed_sites
    mv *removed.txt removed_sites
    
    # And do the same for the log files
    mkdir log_files
    mv *log.txt log_files
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CIAlign: \$( CIAlign --version )
    END_VERSIONS
    """
}
