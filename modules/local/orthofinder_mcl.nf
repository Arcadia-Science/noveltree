process ORTHOFINDER_MCL {
    tag "MCL clustering"

    label 'process_lowcpu'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/orthofinder_2.5.4:0.0.1' :
        '' }"

    // stageInMode = "copy"

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
    # if so, delete the sequence files and other non-essential directories that
    # we will not be using and take up significant, unnecessary space.
    # if [ "$output_directory" == "mcl_test_dataset" ]; then
    #     rm -r OrthoFinder/*/Single_Copy_Orthologue_Sequences/
    #     rm -r OrthoFinder/*/Orthogroup_Sequences/
    #     rm -r OrthoFinder/*/WorkingDirectory/
    #     rm -r OrthoFinder/*/Orthologues/
    # else
    #     dir=\$(pwd)
    #     cd \$(ls -d OrthoFinder/*/WorkingDirectory)
    #     tar -czvf Sequences_ids.tar.gz Sequences_ids
    #     rm -r Sequences_ids
    #     cd \$dir
    # fi

    # Restructure to get rid of the unnecessary "OrthoFinder" directory"
    #mkdir ${output_directory}
    #mv OrthoFinder/Results_Inflation_${mcl_inflation}/ ${output_directory}/
    #rm -r OrthoFinder/
    """
}
