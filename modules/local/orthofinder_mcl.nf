process ORTHOFINDER_MCL {
    tag "$blastdir"
    label 'process_medium'
    container "davidemms/orthofinder"

    input:
    path blastdir
    each val(inflation)

    output:
    path "${fasta}.dmnd", emit: ogs
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def inflation_param = inflation ? "${inflation}" : '1.5'
    """
    orthofinder \\
        -b $blastdir \\
        -I $inflation_param \\
        -M msa -X -os -z \\
        -a $task.cpus \\
        --in $fasta \\
        -d $fasta \\
        $args
        
    mv ${fasta_dir}/OrthoFinder/Results*/WorkingDirectory/* ${fasta_dir}/OrthoFinder/
    rm -r ${fasta_dir}/OrthoFinder/Results*/WorkingDirectory/
    rm -r ${fasta_dir}/OrthoFinder/Results*/
    
    """
}

