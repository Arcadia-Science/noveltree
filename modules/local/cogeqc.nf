process COGEQC {
    tag "Orthogroup Summary"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::bioconductor-uniprot.ws==2.34.0--r41hdfd78af_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-uniprot.ws:2.34.0--r41hdfd78af_0':
        'bioconductor/bioconductor_docker' }"
        
    publishDir(
        path: "${params.outdir}/orthogroup-summaries",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path orthofinder_outdir // Files storing filepaths to where orthofinder results are stored for each inflation parameter
    path prot_annotations // Base filepath to where protein annotations are stored

    output:
    path "*-cogeqc-summary.tsv", emit: og_summary
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    def orthogroups_fpath = "${orthofinder_outdir}"
    //def annotations = "${prot_annotations}"
    """
    # Assess orthogroups inferred using each inflation parameter, summarizing 
    # how well they group proteins with the same domains together, as well as
    # other summary stats like number of ogs with >= 4 species, per-species
    # gene count per-og, etc. 
        
    # store the orthogroup directory as a bash variable using the file storing
    # the filepath used to stall initiation of this module until 
    orthogroups=\$(cat $orthogroups_fpath)
    Rscript $projectDir/bin/cogeqc-summarize-ogs.R \$orthogroups
    
    # Additionally, move the orthogroup inference directories to a new 
    # location, so that these 'test' orthogroups do not interfere with 
    # the full orthogroup inference procedure. 
    #mv ../../../${params.outdir}/orthofinder/Results_Inflation_* ../../../${params.outdir}/orthogroup-summaries/
    #mv ../../../${params.outdir}/orthofinder/MCL*fpath.txt ../../../${params.outdir}/orthogroup-summaries/
    #mv ../../../${params.outdir}/orthofinder/mcl_testing/Results_Inflation_* ../../../${params.outdir}/orthogroup-summaries/
    #mv ../../../${params.outdir}/orthofinder/mcl_testing/MCL*fpath.txt ../../../${params.outdir}/orthogroup-summaries/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cogeqc: \$( cat version.txt | sed "s/\\[1] ‘//g" | sed "s/’//g" )
    END_VERSIONS
    """
}
