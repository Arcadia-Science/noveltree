process ORTHOFINDER_MCL {
    tag "MCL clustering"
    label 'process_lowcpu'
    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder:2.5.4' :
        '' }"

    input:
    each mcl_inflation
    path(blast)
    path(fasta)
    path(db)
    path(sppIDs)
    path(seqIDs)
    val output_directory

    output:
    path("*/Results_Inflation*"), emit: inflation_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    for f in \$(ls TestBlast*)
    do
        mv \$f \$(echo \$f | sed "s/TestBlast/Blast/g")
    done

    orthofinder \\
        -b ./ \\
        -n "Inflation_${mcl_inflation}" \\
        -I $mcl_inflation \\
        -M msa -X -os -z \\
        -a ${task.cpus} \\
        $args

    # Check if we're running an mcl test or not:
    # if so, delete the sequence files, which we will not be using and take up
    # significant, unnecessary space.
    if [ "$output_directory" == "mcl_test_dataset" ]; then
        rm -r OrthoFinder/*/Orthogroup_Sequences/
    fi

    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    mkdir ${output_directory}
    mv OrthoFinder/Results_Inflation_${mcl_inflation}/ ${output_directory}/
    rm -r OrthoFinder/
    """
}
