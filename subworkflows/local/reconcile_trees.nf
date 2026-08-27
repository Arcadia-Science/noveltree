//
// Species tree inference + gene-tree/species-tree reconciliation
//

include { SPECIESRAX          } from '../../modules/local/speciesrax'
include { BUNDLE_SPECIESRAX_INPUTS } from '../../modules/local/bundle_speciesrax_inputs'
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
    reference_chronogram   // optional rooted TimeTree chronogram; empty with explicit outgroups
    speciesrax_selection   // upstream family-selection audit table
    speciesrax_coverage    // upstream per-species selected-family coverage

    main:
    // Join each core family's mapping and tree, then create bounded
    // uncompressed shards. With --strategy SKIP, SpeciesRax never reads an
    // alignment, so excluding alignments avoids staging the largest inputs.
    ch_speciesrax_family_inputs = core_map_links
        .join(core_gene_trees)

    def speciesrax_bundle_size = params.speciesrax_bundle_size as int
    ch_speciesrax_bundle_inputs = ch_speciesrax_family_inputs
        .map { meta, map_link, gene_tree ->
            def matcher = meta.og =~ /^OG(\d+)$/
            if (!matcher.matches()) {
                throw new IllegalArgumentException("Unexpected orthogroup identifier: ${meta.og}")
            }
            def shard_number = (matcher[0][1] as long).intdiv(speciesrax_bundle_size)
            tuple(String.format('%08d', shard_number), meta, map_link, gene_tree)
        }
        .groupTuple()
        .map { shard_id, metas, map_links, gene_trees ->
            def families = (0..<metas.size()).collect { index ->
                [metas[index], map_links[index], gene_trees[index]]
            }
            def ordered = families.sort { left, right -> left[0].og <=> right[0].og }
            tuple(
                shard_id,
                ordered.collect { it[1] },
                ordered.collect { it[2] },
            )
        }

    BUNDLE_SPECIESRAX_INPUTS(ch_speciesrax_bundle_inputs)

    //
    // SPECIESRAX: infer MiniNJ, transfer the trusted root, then estimate final
    // branch lengths/support while preserving that rooted topology.
    //
    SPECIESRAX(
        BUNDLE_SPECIESRAX_INPUTS.out.archive.collect(),
        reference_chronogram,
        species_name_list,
        params.outgroups,
        speciesrax_selection,
        speciesrax_coverage
    )
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
    species_coverage      = GENERAX_PER_SPECIES.out.species_coverage
    generax_nhx           = GENERAX_PER_SPECIES.out.generax_nhx
    labeled_species_tree  = GENERAX_PER_SPECIES.out.labeled_species_tree
    generax_per_spp_gfts  = GENERAX_PER_SPECIES.out.generax_per_spp_gfts
    all_clean_msas        = ch_all_og_clean_msas
}
