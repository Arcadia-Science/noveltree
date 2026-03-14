process GENERAX_PER_SPECIES {
    tag "$meta.og"

    cpus { 16 * task.attempt }
    time { 6.h * task.attempt }
    memory {
        def n = (meta?.n_seq ?: 50) as long
        def L = (meta?.max_len ?: 500) as long
        def clv_bytes = n * L * 640L
        def mpi_overhead = task.cpus * 250L * 1024L * 1024L
        def estimated_gb = Math.max(8L, (long)((clv_bytes * 2.5 + mpi_overhead) / (1024L * 1024L * 1024L)) + 2L)
        def capped_gb = (int) Math.min(estimated_gb, 96L)
        def requested = capped_gb.GB * task.attempt
        def max_mem = params.max_memory as nextflow.util.MemoryUnit
        requested.compareTo(max_mem) > 0 ? max_mem : requested
    }

    stageInMode 'copy' // Must stage in as copy, or OpenMPI will try to contantly read from S3 which causes problems.

    // Exit code 10 = "no valid families" (e.g. invalid starting tree).
    // This is deterministic — retrying won't help. Skip the OG gracefully
    // so one bad gene family doesn't kill the entire pipeline.
    errorStrategy { task.exitStatus == 10 ? 'ignore' : (task.attempt <= 5 ? 'retry' : 'terminate') }
    maxRetries 5

    container 'arcadiascience/generax_56f3ed0:1.1.3'

    storeDir "${params.outdir}/generax/per_species_rates"
    publishDir(
        path: "${params.outdir}/gene_family_trees/reconciled",
        mode: params.publish_dir_mode,
        pattern: "*/*_reconciled_gft.newick",
        saveAs: { fn -> fn.split('/')[-1].replace('_reconciled_gft', '_gps_reconciled') },
    )

    input: // Input is a single large tuple with paths to map-links, tree files, alignments, and the species tree
    tuple val(meta), file(map_link), file(gene_tree), file(alignment), file(species_tree)

    output:
    tuple val(meta), path("${meta.og}/${meta.og}_reconciled_gft.newick")        , emit: generax_per_spp_gfts
    tuple val(meta), path("${meta.og}/${meta.og}_eventCounts.txt")              , emit: event_counts
    tuple val(meta), path("${meta.og}/${meta.og}_speciesEventCounts.txt")       , emit: species_event_counts
    tuple val(meta), path("${meta.og}/${meta.og}_transfers.txt")                , emit: transfer_event_counts
    tuple val(meta), path("${meta.og}/${meta.og}_perSpeciesCoverage.txt")       , emit: species_coverage
    tuple val(meta), path("${meta.og}/${meta.og}_reconciliated.nhx")            , emit: generax_nhx
    path "generax_labeled_species_tree.newick"                                  , emit: labeled_species_tree
    path "${meta.og}/${meta.og}_full_output.tar.gz"                             , emit: archive

    when:
    task.ext.when == null || task.ext.when

    // always gets set as the file itself, excluding the path
    script:
    def args = task.ext.args ?: ''
    def og   = "${meta.og}"
    """
    # Recode selenocysteine as a gap character:
    # RAxML-NG (used under the hood by SpeciesRax and
    # GeneRax) cannot handle these. Even if rare,
    # their inclusion leads a number of gene families
    # to be excluded from analyses.
    sed -E -i '/>/!s/U/X/g' *.fa

    # Do the same for Pyrrolysine
    sed -E -i '/>/!s/O/X/g' *.fa

    # Populate the family file for this gene family for the
    # analysis with GeneRax
    # We will be using LG+G4+F for all gene families
    echo "[FAMILIES]" > ${og}.family
    echo "- ${og}" >> ${og}.family
    echo "starting_gene_tree = ${gene_tree}" >> ${og}.family
    echo "mapping = ${map_link}" >> ${og}.family
    echo "alignment = $alignment" >> ${og}.family
    echo "subst_model = LG+G4+F" >> ${og}.family

    mpiexec \\
        -np ${task.cpus} \\
        --allow-run-as-root \\
        --use-hwthread-cpus \\
        generax \\
        --species-tree $species_tree \\
        --families ${og}.family \\
        --per-species-rates \\
        --prefix $og \\
        --reconciliation-samples 100 \\
        $args

    # Clean up
    rm -fr $og/gene_optimization_*

    # Use the events tree as the reconciled gene tree — it's the same topology and branch
    # lengths as geneTree.newick but includes S/D/T internal node labels from reconciliation
    cp "$og/reconciliations/${og}_events.newick" $og/results/$og/${og}_reconciled_gft.newick

    # Rename perSpeciesCoverage.txt to include orthogroup prefix (prevents file name collisions downstream)
    mv "$og/perSpeciesCoverage.txt" "$og/${og}_perSpeciesCoverage.txt"

    # And move the reconciliation transfer samples into a subdirectory, archive, and compress.
    mkdir $og/reconciliations/reconciliation_transfer_samples/
    mv $og/reconciliations/*_*_transfers.txt $og/reconciliations/reconciliation_transfer_samples/
    tar -czvf $og/reconciliations/reconciliation_transfer_samples.tar.gz $og/reconciliations/reconciliation_transfer_samples/
    rm -r $og/reconciliations/reconciliation_transfer_samples/

    # Extract key files to working directory
    cp $og/results/$og/${og}_reconciled_gft.newick .
    cp $og/reconciliations/${og}_eventCounts.txt .
    cp $og/reconciliations/${og}_speciesEventCounts.txt .
    cp $og/reconciliations/${og}_transfers.txt .
    cp $og/reconciliations/${og}_reconciliated.nhx .
    mv $og/${og}_perSpeciesCoverage.txt .

    # Extract GeneRax-labeled species tree (has Node_X_Y_0 internal labels, identical across all OGs)
    cp $og/species_trees/inferred_species_tree.newick generax_labeled_species_tree.newick

    # Archive full GeneRax output, then replace with flat structure
    tar -czf ${og}_full_output.tar.gz $og/
    rm -rf $og/
    mkdir $og
    mv ${og}_reconciled_gft.newick $og/
    mv ${og}_eventCounts.txt $og/
    mv ${og}_speciesEventCounts.txt $og/
    mv ${og}_transfers.txt $og/
    mv ${og}_perSpeciesCoverage.txt $og/
    mv ${og}_reconciliated.nhx $og/
    mv ${og}_full_output.tar.gz $og/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        generax: \$(generax --version | head -n1 | sed 's/.*GeneRax //')
    END_VERSIONS
    """
}
