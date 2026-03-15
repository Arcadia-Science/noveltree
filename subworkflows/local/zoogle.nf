//
// Zoogle: physicochemical properties, time calibration, dating, and phylo-distance analysis
//

include { PROTEIN_PROPERTIES          } from '../../modules/local/protein_properties'
include { BUILD_REFERENCE_CHRONOGRAM  } from '../../modules/local/build_reference_chronogram'
include { TIME_CALIBRATE_SPECIES_TREE } from '../../modules/local/time_calibrate_species_tree'
include { DATE_GENE_FAMILY_TREES      } from '../../modules/local/date_gene_family_trees'
include { ZOOGLE_ANALYSIS             } from '../../modules/local/zoogle_analysis'

workflow ZOOGLE {
    take:
    spptree_fas            // [ val(meta), path(fasta) ]
    genetree_fas           // [ val(meta), path(fasta) ]
    all_clean_msas         // [ val(meta), path(msa) ]
    generax_per_spp_gfts   // [ val(meta), path(tree) ]
    speciesrax_tree        // path
    species_name_list      // value channel
    ref_species            // string value

    main:
    //
    // PROTEIN_PROPERTIES: physicochemical properties and AA composition for all gene families
    //
    ch_all_og_original_fas = spptree_fas.concat(genetree_fas)
    ch_physchem_input = ch_all_og_original_fas.join(all_clean_msas)
    PROTEIN_PROPERTIES(ch_physchem_input)

    //
    // Time calibration: either auto-build from TimeTree.org or use user-provided tree
    //
    if (!params.reference_time_tree || params.reference_time_tree == 'none') {
        ch_species_names_file = species_name_list
            .collectFile(name: 'species_names.txt', newLine: true)
        BUILD_REFERENCE_CHRONOGRAM(ch_species_names_file, params.ncbi_email)
        ch_reference_tree = BUILD_REFERENCE_CHRONOGRAM.out.chronogram
    } else {
        ch_reference_tree = file(params.reference_time_tree)
    }

    TIME_CALIBRATE_SPECIES_TREE(
        speciesrax_tree,
        ch_reference_tree,
        params.time_calibration_method,
        params.age_bracket
    )

    //
    // DATE_GENE_FAMILY_TREES: time-calibrate gene family trees using
    // reconciliation-filtered (speciation-only) calibrations
    //
    ch_dating_input = generax_per_spp_gfts
        .join(all_clean_msas)
        .combine(TIME_CALIBRATE_SPECIES_TREE.out.calibrated_tree)

    DATE_GENE_FAMILY_TREES(ch_dating_input, params.max_treepl_tips)

    //
    // ZOOGLE_ANALYSIS: phylogenetic distance analysis with validation filter
    //
    ch_zoogle_input = DATE_GENE_FAMILY_TREES.out.dated_gft
        .join(PROTEIN_PROPERTIES.out.summary_stats)
        .filter { meta, tree, props_file ->
            // Validate gene family has sufficient proteins for phylo-dist analysis
            // Read CSV and extract protein IDs (first column, skip header)
            def lines = props_file.readLines()
            def proteinIds = lines.drop(1).collect { it.split(',')[0] }

            if (proteinIds.size() == 0) {
                log.info "Skipping ${meta.og}: No proteins found in CSV"
                return false
            }

            // Count reference species proteins
            def refCount = proteinIds.count { it.startsWith("${ref_species}_") }

            // Count non-reference proteins
            def nonrefCount = proteinIds.size() - refCount

            // Count proteins per non-reference species
            // NOTE: Must match R script's species extraction logic (line 145 of protein_distance_calculation_functions.R)
            // Protein labels have format: Species_name_ProteinID
            // Species names are extracted by removing the last underscore-delimited segment
            def nonrefProteinsBySpecies = proteinIds
                .findAll { !it.startsWith("${ref_species}_") }
                .collect { it.replaceFirst(/_[^_]+$/, '') }  // Extract species name (remove protein ID after last underscore)
                .countBy { it }  // Map of species -> count

            // Count unique non-reference species
            def nonrefSpeciesCount = nonrefProteinsBySpecies.size()

            // Count how many non-reference species have at least 2 proteins
            // (Wilcoxon test requires at least 2 observations per group)
            def speciesWithEnoughProteins = nonrefProteinsBySpecies.count { species, count -> count >= 2 }

            // Apply validation criteria (maps directly to the 3 observed errors)
            def isValid = (refCount >= 1) && (nonrefCount >= 2) && (speciesWithEnoughProteins >= 2)

            return isValid
        }

    ZOOGLE_ANALYSIS(ch_zoogle_input, ref_species)

    emit:
    calibrated_tree = TIME_CALIBRATE_SPECIES_TREE.out.calibrated_tree
}
