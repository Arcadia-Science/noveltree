process WITCH {
    tag "$meta.og"
    label 'process_witch'

    container "${ workflow.containerEngine == 'docker' ? 'arcadiascience/witch_0.3.0:1.0.0' :
        '' }"
    // TODO: address this issue (permission related errors) in future release
    containerOptions = "--user root"

    stageInMode = 'copy'
    publishDir(
        path: "${params.outdir}/witch_alignments",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("**_witch.fa")         , emit: msas
    tuple val(meta), path("**_witch_cleaned.fa") , emit: cleaned_msas, optional: true
    tuple val(meta), path("**_map.link")         , emit: map_link, optional: true
    path("*")                                    , emit: results
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def og      = "${meta.og}"
    def min_len = params.min_ungapped_length ?: '20'
    """
    # Scale down the number of threads allocated to MAFFT based on the number of
    # task attempts (i.e. reduce the memory usage of MAFFT if MAGUS keeps 
    # running out of memory for this reason). Will slow down MAGUS somewhat, 
    # but will increase the chance of the process completing successfully. 
    mafft_threads=\$(( 9 - ${task.attempt} ))
    sed -i "s/numthreadsit=8/numthreadsit=\$mafft_threads/g" /WITCH/*/*/*/*/*/bin/mafft 
    sed -i "s/numthreads -lt 8/numthreads -lt \$mafft_threads/g" /WITCH/*/*/*/*/*/bin/mafft 
    
    # If we are resuming a run, do some cleanup:
    if [ -d "alignments/" ]; then
        rm -rf alignments/
    fi

    # Be sure to remove any non-standard amino acid codes in the input sequences, as this
    # can cause errors downstream and in parsing.
    sed -E -i '/>/!s/U//g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O//g' ${fasta} # pyrrolysine

    python3 /WITCH/witch.py \\
        -i ${fasta} \\
        -d alignments \\
        -t ${task.cpus} \\
        --graphtraceoptimize true \\
        --molecule amino \\
        $args

    # By default, witch only masks (and removes) singleton columns - we still would
    # like to make sure that we are removing sequences with too many positions removed.
    # So, here we're using awk to remove sequences with fewer than params.min_ungapped_length
    # AA remaining once masked.
    awk -v N=${min_len} -F "" \
        'BEGIN { getline; header=\$0; seq="" } \
        !/^>/ { for (i=1; i<=NF; i++) if (\$i != "-") s++ } \
        /^>/ { if (s >= N || seq == "") { if (header != "") print header; if (seq != "") print seq } header=\$0; seq=""; s=0 } \
        !/^>/ { seq = seq \$0 } \
        END { if (s >= N) { print header; print seq } }' \
        alignments/merged.fasta.masked > tmp.fasta

    # And remove any columns that are now comprised exclusively of gaps following the exclusion
    # of (if any) sequences in the above step.
    awk 'BEGIN {seq_count=0} \
        /^>/ {seq_count++; headers[seq_count]=\$0; next} \
        {sequences[seq_count]=sequences[seq_count]\$0} \
        END {for(i=1;i<=length(sequences[1]);i++){ \
            column=""; \
            for(j=1;j<=seq_count;j++){column=column substr(sequences[j],i,1)} \
            if(column!~/^-+\$/){ \
              for(j=1;j<=seq_count;j++){new_sequences[j]=new_sequences[j] substr(sequences[j],i,1)}}\
          } \
          for(i=1;i<=seq_count;i++){print headers[i]; print new_sequences[i]} \
        }' tmp.fasta > final_masked.fasta

    # Reorganize results for publishing
    mkdir original_alignments
    mkdir cleaned_alignments
    mv alignments/merged.fasta original_alignments/${og}_witch.fa
    mv final_masked.fasta cleaned_alignments/${og}_witch_cleaned.fa
    rm -r alignments/ && rm tmp.fasta

    # In the rare case that this filtering reduces sequences down to < 4
    # sequences, delete the output cleaned alignments to exclude them from
    # downstream phylogenetic analyses.
    n_remain=\$(grep ">" cleaned_alignments/${og}_witch_cleaned.fa | wc -l)
    if [ \$n_remain -lt 4 ]; then
        rm cleaned_alignments/${og}_witch_cleaned.fa
    else
        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        mkdir species_protein_maps
        grep ">" cleaned_alignments/${og}_witch_cleaned.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > species_protein_maps/${og}_map.link
        rm prot && rm spp
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        witch: v\$(python3 /WITCH/witch.py -v | cut -f2 -d " ")
    END_VERSIONS
    """
}
