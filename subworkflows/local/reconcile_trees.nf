//
// Species tree inference + gene-tree/species-tree reconciliation
//

if (params.outgroups != 'none') {
    include { ASTEROID } from '../../modules/local/asteroid'
}
include { SPECIESRAX          } from '../../modules/local/speciesrax'
include { GENERAX_PER_SPECIES } from '../../modules/local/generax_per_species'
if (params.generax_per_family) {
    include { GENERAX_PER_FAMILY } from '../../modules/local/generax_per_family'
}

workflow RECONCILE_TREES {
    take:
    species_name_list      // value channel: list of species names
    core_gene_trees        // [ val(meta), path(tree) ]
    core_map_links         // [ val(meta), path(map_link) ]
    core_clean_msas        // [ val(meta), path(msa) ]
    remaining_gene_trees   // [ val(meta), path(tree) ]
    remaining_map_links    // [ val(meta), path(map_link) ]
    remaining_clean_msas   // [ val(meta), path(msa) ]

    main:
    // Collect core-set lists for SPECIESRAX
    core_gene_tree_list    = core_gene_trees.collect { it[1] }
    core_og_maplink_list   = core_map_links.collect { it[1] }
    core_og_clean_msa_list = core_clean_msas.collect { it[1] }

    //
    // ASTEROID: infer initial unrooted species tree (optional, when outgroups provided)
    //
    if (params.outgroups != "none") {
        ASTEROID(species_name_list, core_gene_tree_list, params.outgroups)
            .rooted_spp_tree
            .set { ch_asteroid }
    } else {
        ch_asteroid = Channel.value("none")
    }

    //
    // SPECIESRAX: infer rooted species tree with gene-tree/species-tree reconciliation
    //
    SPECIESRAX(core_og_maplink_list, core_gene_tree_list, core_og_clean_msa_list, ch_asteroid)
        .speciesrax_tree
        .set { ch_speciesrax }

    // Merge core + remaining for GeneRax
    ch_all_map_links     = core_map_links.concat(remaining_map_links)
    ch_all_gene_trees    = core_gene_trees.concat(remaining_gene_trees)
    ch_all_og_clean_msas = core_clean_msas.concat(remaining_clean_msas)

    // Join so each gene family is processed asynchronously, combined with species tree
    ch_generax_input = ch_all_map_links
        .join(ch_all_gene_trees)
        .join(ch_all_og_clean_msas)
        .combine(ch_speciesrax)

    // GENERAX_PER_FAMILY: full mode only
    if (params.generax_per_family) {
        GENERAX_PER_FAMILY(ch_generax_input)
    }

    // GENERAX_PER_SPECIES: both modes
    GENERAX_PER_SPECIES(ch_generax_input)

    emit:
    speciesrax_tree       = ch_speciesrax
    event_counts          = GENERAX_PER_SPECIES.out.event_counts
    species_event_counts  = GENERAX_PER_SPECIES.out.species_event_counts
    transfer_event_counts = GENERAX_PER_SPECIES.out.transfer_event_counts
    species_coverage      = GENERAX_PER_SPECIES.out.species_coverage
    generax_nhx           = GENERAX_PER_SPECIES.out.generax_nhx
    labeled_species_tree  = GENERAX_PER_SPECIES.out.labeled_species_tree
    generax_per_spp_gfts  = GENERAX_PER_SPECIES.out.generax_per_spp_gfts
    all_clean_msas        = ch_all_og_clean_msas
}
