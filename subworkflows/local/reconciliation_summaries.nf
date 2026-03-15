//
// Phylogenetic profiles, HGT summaries, and hierarchical ortholog group parsing
//

include { PHYLO_PROFILES       } from '../../modules/local/phylo_profiles'
include { MERGE_PHYLO_PROFILES } from '../../modules/local/merge_phylo_profiles'
include { PARSE_PHYLOHOGS      } from '../../modules/local/parse_phylohogs'

workflow RECONCILIATION_SUMMARIES {
    take:
    event_counts           // [ val(meta), path(tsv) ]
    species_event_counts   // [ val(meta), path(tsv) ]
    transfer_event_counts  // [ val(meta), path(tsv) ]
    species_coverage       // [ val(meta), path(tsv) ]
    inflation_dir          // path
    generax_nhx            // [ val(meta), path(nhx) ]
    labeled_species_tree   // path

    main:
    //
    // PHYLO_PROFILES: generate per-gene-family phylogenetic profiles
    //
    ch_phylo_profiles_input = event_counts
        .join(species_event_counts)
        .join(transfer_event_counts)
        .join(species_coverage)

    PHYLO_PROFILES(ch_phylo_profiles_input, inflation_dir)

    // Merge per-OG phylo profile outputs.
    // collectFile() concatenates TSVs natively in Nextflow, avoiding the need
    // to stage hundreds/thousands of small files from S3 into a single process.
    PHYLO_PROFILES.out.duplication_count
        .collectFile(name: 'duplication_count_per_species_per_gene_family.tsv',
                     storeDir: "${params.outdir}/gene_family_evolution", keepHeader: true)
    PHYLO_PROFILES.out.loss_count
        .collectFile(name: 'loss_count_per_per_species_gene_family.tsv',
                     storeDir: "${params.outdir}/gene_family_evolution", keepHeader: true)
    PHYLO_PROFILES.out.speciation_count
        .collectFile(name: 'speciation_count_per_species_per_gene_family.tsv',
                     storeDir: "${params.outdir}/gene_family_evolution", keepHeader: true)
    PHYLO_PROFILES.out.transfer_donor_count
        .collectFile(name: 'transfer_donor_count_per_species_per_gene_family.tsv',
                     storeDir: "${params.outdir}/gene_family_evolution", keepHeader: true)
    PHYLO_PROFILES.out.transfer_recipient_count
        .collectFile(name: 'transfer_recipient_count_per_species_per_gene_family.tsv',
                     storeDir: "${params.outdir}/gene_family_evolution", keepHeader: true)

    // HGT counts are emitted in long format (donor, recipient, count).
    // collectFile() concatenates them, then a lightweight process pivots
    // the single file into the final species × species matrix.
    ch_hgt_long = PHYLO_PROFILES.out.hgt_counts_long
        .collectFile(name: 'hgt_counts_long_all.tsv', keepHeader: true)
    MERGE_PHYLO_PROFILES(ch_hgt_long)

    //
    // PARSE_PHYLOHOGS: extract ortholog/paralog/xenolog pairs and HOG membership
    //
    PARSE_PHYLOHOGS(generax_nhx, labeled_species_tree)

    emit:
    hgt_summed_count = MERGE_PHYLO_PROFILES.out.hgt_summed_count
}
