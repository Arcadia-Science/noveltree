process DIAMOND_BLASTP {
    // Modified from nf-core to flexibly handle application to both complete and mcl-test datasets

    tag "$meta.id"
    label 'process_diamond'

    conda (params.enable_conda ? "bioconda::diamond=2.0.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.0.15--hb97b32f_0' :
        'quay.io/biocontainers/diamond:2.0.15--hb97b32f_0' }"

    input:
    tuple val(meta), path("*")
    file(fasta)
    each file(db)
    val output_extension
    val mcl_test

    output:
    path('*.blast*')    , optional: true, emit: blast
    path('*.xml*')      , optional: true, emit: xml
    path('*.txt*')      , optional: true, emit: txt
    path('*.daa*')      , optional: true, emit: daa
    path('*.sam*')      , optional: true, emit: sam
    path('*.tsv*')      , optional: true, emit: tsv
    path('*.paf*')      , optional: true, emit: paf
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def testing_mcl = mcl_test.equals('true') ? "${mcl_test}" : "false"
    switch ( output_extension ) {
        case "blast": outfmt = 0; break
        case "xml": outfmt = 5; break
        case "txt": outfmt = 6; break
        case "daa": outfmt = 100; break
        case "sam": outfmt = 101; break
        case "tsv": outfmt = 102; break
        case "paf": outfmt = 103; break
        default:
            outfmt = '6';
            output_extension = 'txt';
            log.warn("Unknown output file format provided (${output_extension}): selecting DIAMOND default of tabular BLAST output (txt)");
            break
    }
    """
    # Get the species name against which we're querying
    spp_query=\$(echo $fasta | sed "s/Species//g" | sed 's/.fa//g' | sed 's|.*/||g')
    sbb_db=\$(echo $db | sed "s/diamondDBSpecies//g" | sed 's/.dmnd//g' | sed 's|.*/||g')

    if [ "$testing_mcl" == "true" ]; then
        out_name="TestBlast\${spp_query}_\${sbb_db}.${output_extension}"
    else
        out_name="Blast\${spp_query}_\${sbb_db}.${output_extension}"
    fi

    diamond \\
        blastp \\
        --out \$out_name \\
        --outfmt ${outfmt} \\
        --threads ${task.cpus} \\
        --query $fasta \\
        --compress 1 \\
        --db $db \\
        $args

    cat <<-END_VERSIONS >> versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
