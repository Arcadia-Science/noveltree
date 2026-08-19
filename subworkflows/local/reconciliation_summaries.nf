//
// Phylogenetic profiles and hierarchical ortholog group parsing
//

include { PHYLO_PROFILES       } from '../../modules/local/phylo_profiles'
include { PARSE_PHYLOHOGS      } from '../../modules/local/parse_phylohogs'

workflow RECONCILIATION_SUMMARIES {
    take:
    species_event_counts   // [ val(meta), path(tsv) ]
    species_coverage       // [ val(meta), path(tsv) ]
    generax_nhx            // [ val(meta), path(nhx) ]
    labeled_species_tree   // path

    main:
    //
    // PHYLO_PROFILES: generate per-gene-family phylogenetic profiles
    //
    ch_phylo_profiles_input = species_event_counts.join(species_coverage)

    PHYLO_PROFILES(ch_phylo_profiles_input)

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

    //
    // PARSE_PHYLOHOGS: extract ortholog/paralog pairs and HOG membership
    //
    PARSE_PHYLOHOGS(generax_nhx, labeled_species_tree.first())

    emit:
    orthologs        = PARSE_PHYLOHOGS.out.orthologs
    paralogs         = PARSE_PHYLOHOGS.out.paralogs
}
