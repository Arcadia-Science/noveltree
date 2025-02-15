process ORTHOFINDER_PREP {
    tag "Prepping data for OrthoFinder"
    label 'process_low'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder_2.5.4:1.0.0' :
        '' }"

    publishDir(
        path: "${params.outdir}/orthofinder",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file(fasta)
    val output_directory

    output:
    path "**.dmnd"           , emit: diamonds
    path "**.fa"             , emit: fastas
    path "**SequenceIDs.txt" , emit: seqIDs
    path "**SpeciesIDs.txt"  , emit: sppIDs
    path "versions.yml"      , emit: versions

    script:
    def sppid_protid_delim = "${params.sppid_protid_delim}"
    """
    # Make sure the delimiter is actually present (once) in the sequence ID
    # Use a representative fasta file to find out
    rep_fa=\$(ls ./*.fa* | head -n1)
    seqid=\$(head -n1 \$rep_fa | cut -f1 -d" ")
    in_seqid=\$(echo \$seqid | awk -F'$sppid_protid_delim' '{print NF-1}')
    
    # If the delimiter is not actually in the seqid, stop the workflow.
    if [ \$in_seqid == 0 ]; then
        echo "ERROR: The species and protein ID delimiter is not present in the sequence IDs!"
        echo "Double check that you are using the correct delimiter for your data as specified using the 'sppid_protid_delim' workflow parameter. This delimiter must only occur once, separating the species ID and protein ID" 
        echo "You are using the the following delimiter: '$sppid_protid_delim'"
        exit 1
    else
        # The fasta directory depends on whether we're running the mcl testing or not.
        orthofinder \\
            -f ./ \\
            -t ${task.cpus} \\
            -op > tmp
            
        mkdir ${output_directory} && mv OrthoFinder/ ${output_directory}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder --versions | head -n2 | tail -n1 | sed "s/OrthoFinder version //g" | sed "s/ Copyright (C) 2014 David Emms//g")
    END_VERSIONS
    """
}
