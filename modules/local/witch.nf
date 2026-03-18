process WITCH {
    tag "$meta.og"

    cpus { Math.min( 16, params.max_cpus as int ) }
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 1000) as long
        def L = (meta?.max_len ?: 500) as long
        // WITCH internally runs MAGUS which runs MAFFT on a backbone of up to 1000 seqs.
        // MAFFT E-INS-i is O(N^2 * L) — the dominant memory cost.
        // Then WITCH runs HMMER searches (O(N * L)) and stores extended alignment (O(N * L)).
        def backbone_n = Math.min(n, 1000L)
        def mafft_gb = (long)(backbone_n * backbone_n * L * 8L / (1024L * 1024L * 1024L))
        def hmm_gb = (long)(n * L * 16L / (1024L * 1024L * 1024L))
        // Base overhead: Python, MAGUS, FastTree, MCL, HMMER processes
        def estimated_gb = Math.max(32L, mafft_gb + hmm_gb + 8L)
        def capped_gb = (int) Math.min(estimated_gb, 256L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    container 'arcadiascience/witch_1.0.10:1.0.0'
    // WITCH may write to its install directory, requiring writable container filesystem
    containerOptions = workflow.containerEngine == 'docker' ? '--user root' : \
        (workflow.containerEngine == 'singularity' ? '--writable-tmpfs' : '')

    stageInMode = 'copy'
    storeDir "${params.outdir}/alignments"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("trimmed/${fasta.baseName}_witch_cleaned.fa"), emit: msas, optional: true
    tuple val(meta), path("trimmed/species_protein_maps/${fasta.baseName}_map.link"), emit: map_link, optional: true
    path "original/${fasta.baseName}_witch.fa"   , emit: original_alignments

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def og      = "${meta.og}"
    def min_len = params.min_ungapped_length ?: '20'
    def min_seq = params.min_num_seq_per_og
    def min_spp = params.min_num_spp_per_og
    """
    # If we are resuming a run, do some cleanup:
    if [ -d "alignments/" ]; then
        rm -rf alignments/
    fi

    # Be sure to remove any non-standard amino acid codes in the input sequences, as this
    # can cause errors downstream and in parsing.
    sed -E -i '/>/!s/U/X/g' ${fasta} # selenocysteine
    sed -E -i '/>/!s/O/X/g' ${fasta} # pyrrolysine

    witch-msa \\
        -i ${fasta} \\
        -d alignments \\
        -t ${task.cpus} \\
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
        alignments/aligned.masked.fasta > tmp.fasta

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

    # Reorganize results for storeDir
    mkdir -p original
    mkdir -p trimmed
    mv alignments/aligned.fasta original/${og}_witch.fa
    mv final_masked.fasta trimmed/${og}_witch_cleaned.fa
    rm -r alignments/ && rm tmp.fasta

    # Verify the trimmed alignment still meets minimum sequence/species thresholds.
    n_seq=\$(grep -c ">" trimmed/${og}_witch_cleaned.fa || true)
    n_spp=\$(grep ">" trimmed/${og}_witch_cleaned.fa | sed "s/>//" | sed "s/_[^_]*\$//" | sort -u | wc -l | tr -d ' ')
    if [ "\$n_seq" -lt "$min_seq" ] || [ "\$n_spp" -lt "$min_spp" ]; then
        rm trimmed/${og}_witch_cleaned.fa
    else
        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        mkdir -p trimmed/species_protein_maps
        grep ">" trimmed/${og}_witch_cleaned.fa | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > trimmed/species_protein_maps/${og}_map.link
        rm prot && rm spp
    fi
    """
}
