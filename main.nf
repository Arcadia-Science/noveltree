#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Arcadia-Science/noveltree
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Arcadia-Science/noveltree
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate input parameters
WorkflowMain.initialize(workflow, params, log)
def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [params.input]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input samplesheet not specified!'
}
if (params.mcl_inflation) {
    mcl_inflation = params.mcl_inflation.toString().split(",").collect { it.trim() }
} else {
    exit 1, 'MCL Inflation parameter(s) not specified!'
}
// Check if zoogle mode has required parameters for time calibration
if (params.zoogle && (!params.reference_time_tree || params.reference_time_tree == 'none')) {
    if (!params.ncbi_email || params.ncbi_email == 'none') {
        exit 1, 'Zoogle mode without --reference_time_tree requires --ncbi_email to auto-build a reference chronogram from TimeTree.org. Please provide --ncbi_email or --reference_time_tree.'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW
//
include { INPUT_CHECK                           } from './subworkflows/local/input_check'
include { MCL_INFLATION_SELECTION               } from './subworkflows/local/mcl_inflation_selection'
include { INFER_TREES as INFER_SPECIES_TREES    } from './subworkflows/local/infer_trees'
include { INFER_TREES as INFER_REMAINING_TREES  } from './subworkflows/local/infer_trees'

//
// MODULE
//
// Modules being run twice (for MCL testing and full analysis)
// needs to be included twice under different names.
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_ALL  } from './modules/local/orthofinder_prep'
include { FILTER_ORTHOGROUPS                        } from './modules/local/filter_orthogroups'
include { ASTEROID                                  } from './modules/local/asteroid'
include { SPECIESRAX                                } from './modules/local/speciesrax'
include { TIME_CALIBRATE_SPECIES_TREE               } from './modules/local/time_calibrate_species_tree'
include { GENERAX_PER_SPECIES                       } from './modules/local/generax_per_species'
include { DATE_GENE_FAMILY_TREES                    } from './modules/local/date_gene_family_trees'
include { ORTHOFINDER_PHYLOHOGS                     } from './modules/local/orthofinder_phylohogs'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_ALL    } from './modules/local/orthofinder_mcl'
include { PHYLO_PROFILES                            } from './modules/local/phylo_profiles'
include { MERGE_PHYLO_PROFILES                      } from './modules/local/merge_phylo_profiles'
include { PHYSICOCHEMICAL_PROPS                     } from './modules/local/physicochemical_props'
include { PHYLO_DIST                                } from './modules/local/phylo_dist'
include { BUILD_REFERENCE_CHRONOGRAM                } from './modules/local/build_reference_chronogram'

if (params.generax_per_family) {
    include { GENERAX_PER_FAMILY                    } from './modules/local/generax_per_family'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE-MODIFIED MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE
//
include { DIAMOND_BLASTP as DIAMOND_BLASTP_ALL      } from './modules/nf-core-modified/diamond_blastp'

if (params.busco) {
    include { BUSCO as BUSCO_SHALLOW                } from './modules/nf-core-modified/busco'
    include { BUSCO as BUSCO_BROAD                  } from './modules/nf-core-modified/busco'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define a function to instantiate a meta.map to correspond gene family names
// with all outputs we'll be producing
// Function to get list of [meta, [file]]
def create_og_channel(Object inputs) {
    // If the input is a string, convert it to a list containing a single element
    if (inputs instanceof String) {
        inputs = [inputs]
    }
    // create list of maps
    def metaList = []
    inputs.each { input ->
        def meta = [:]
        meta.og = input
        metaList.add(meta)
    }
    return metaList
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// WORKFLOW: Run main Arcadia-Science/noveltree analysis pipeline
//
workflow NOVELTREE {
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_all_data = INPUT_CHECK(ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    species_name_list = ch_all_data.complete_prots.collect { it[0].id }
    complete_prots_list = ch_all_data.complete_prots.collect { it[1] }

    //
    // Running steps to find the best mcl_inflation parameter value.
    // These steps will only run if more than one value was provided.
    //
    if (mcl_inflation.size() > 1) {
        // Use MCL_INFLATION_SELECTION subworkflow for both modes
        MCL_INFLATION_SELECTION(
            ch_all_data.mcl_test_prots,
            ch_all_data.annotation_prots,
            mcl_inflation
        )
        ch_best_inflation = MCL_INFLATION_SELECTION.out.best_inflation
        ch_versions = ch_versions.mix(MCL_INFLATION_SELECTION.out.versions)
    } else {
        ch_best_inflation = Channel.of(mcl_inflation.first())
    }

    //
    // MODULE: Run BUSCO (full mode only)
    // Split up into shallow and broad scale runs, since downstream modules
    // do not use these outputs, so multiple busco runs may be conducted
    // simultaneously
    //
    if (params.busco) {
        // Shallow taxonomic scale:
        BUSCO_SHALLOW(
            ch_all_data.complete_prots.filter{ it[0].shallow_db != "NA" },
            "shallow",
            [],
            []
        )

        // Broad taxonomic scale (Eukaryotes)
        BUSCO_BROAD(
            ch_all_data.complete_prots.filter{ it[0].broad_db != "NA" },
            "broad",
            [],
            []
        )
    }

    //
    // MODULE: Prepare directory structure and fasta files according to
    //         OrthoFinder's preferred format for downstream MCL clustering
    //
    ORTHOFINDER_PREP_ALL(complete_prots_list, "complete_dataset")
    ch_versions = ch_versions.mix(ORTHOFINDER_PREP_ALL.out.versions)

    //
    // MODULE: All-v-All diamond/blastp
    //
    // For the full dataset, to be clustered into orthogroups using
    // the best inflation parameter.
    DIAMOND_BLASTP_ALL(
        ch_all_data.complete_prots,
        ORTHOFINDER_PREP_ALL.out.fastas.flatten(),
        ORTHOFINDER_PREP_ALL.out.diamonds.flatten(),
        "txt",
        "false"
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP_ALL.out.versions)

    // Using this best-performing inflation parameter, infer orthogroups for
    // all samples.
    ORTHOFINDER_MCL_ALL(
        ch_best_inflation,
        DIAMOND_BLASTP_ALL.out.txt.collect(),
        ORTHOFINDER_PREP_ALL.out.fastas,
        ORTHOFINDER_PREP_ALL.out.diamonds,
        ORTHOFINDER_PREP_ALL.out.sppIDs,
        ORTHOFINDER_PREP_ALL.out.seqIDs,
        "complete_dataset"
    )

    //
    // MODULE: FILTER_ORTHOGROUPS
    // Subset orthogroups based on their copy number and distribution
    // across species and taxonomic group.
    // The conservative subset will be used for species tree inference,
    // and the remainder will be used to infer gene family trees only.
    FILTER_ORTHOGROUPS(
        INPUT_CHECK.out.complete_samplesheet,
        ORTHOFINDER_MCL_ALL.out.inflation_dir,
        params.min_num_seq_per_og,
        params.min_num_spp_per_og,
        params.min_prop_spp_for_spptree,
        params.max_copy_num_spp_tree
    )

    // Create meta maps for the two sets by just providing the simple name of each orthogroup:
    spptree_og_names = FILTER_ORTHOGROUPS.out.spptree_fas.map { file -> file.simpleName }
    spptree_og_map = spptree_og_names.map { create_og_channel(it) }.flatten()
    genetree_og_names = FILTER_ORTHOGROUPS.out.genetree_fas.map { file -> file.simpleName }
    genetree_og_map = genetree_og_names.map { create_og_channel(it) }.flatten()

    // And now create the tuple of these output fastas paired with the meta map
    ch_spptree_fas = spptree_og_map.merge(FILTER_ORTHOGROUPS.out.spptree_fas.flatten())
    ch_genetree_fas = genetree_og_map.merge(FILTER_ORTHOGROUPS.out.genetree_fas.flatten())

    //
    // TREE INFERENCE: Alignment → Trimming → Phylogeny
    // Different approaches for full vs simplified modes
    //

    INFER_SPECIES_TREES(ch_spptree_fas)
    ch_versions = ch_versions.mix(INFER_SPECIES_TREES.out.versions)

    INFER_REMAINING_TREES(ch_genetree_fas)
    ch_versions = ch_versions.mix(INFER_REMAINING_TREES.out.versions)

    // Set output channels for downstream use
    ch_core_og_maplinks = INFER_SPECIES_TREES.out.map_link
    ch_rem_og_maplinks = INFER_REMAINING_TREES.out.map_link
    ch_core_og_clean_msas = INFER_SPECIES_TREES.out.cleaned_msas
    ch_rem_og_clean_msas = INFER_REMAINING_TREES.out.cleaned_msas
    ch_core_gene_trees = INFER_SPECIES_TREES.out.phylogeny
    ch_rem_gene_trees = INFER_REMAINING_TREES.out.phylogeny

    core_og_maplink_list = ch_core_og_maplinks.collect { it[1] }
    core_og_clean_msa_list = ch_core_og_clean_msas.collect { it[1] }

    // Create a channel/list (no tuple) of just the core trees used by Asteroid
    core_gene_tree_list = ch_core_gene_trees.collect { it[1] }

    // The following two steps will just be done for the core set of
    // orthogroups that will be used to infer the species tree
    //
    // MODULE: ASTEROID
    // Alrighty, now let's infer an intial, unrooted species tree using Asteroid
    //
    ASTEROID(species_name_list, core_gene_tree_list, params.outgroups)
        .rooted_spp_tree
        .set { ch_asteroid }
    ch_versions = ch_versions.mix(ASTEROID.out.versions)

    // If no outgroups are provided (and thus no rooted species tree output
    // by Asteroid), define ch_asteroid as a null/empty channel
    if (params.outgroups == "none") {
        ch_asteroid = Channel.value("none")
    }

    //
    // MODULE: SPECIESRAX
    // Now infer the rooted species tree with SpeciesRax,
    // reconcile gene family trees, and infer per-family
    // rates of gene-family duplication, transfer, and loss
    //
    SPECIESRAX(core_og_maplink_list, core_gene_tree_list, core_og_clean_msa_list, ch_asteroid)
        .speciesrax_tree
        .set { ch_speciesrax }
    ch_versions = ch_versions.mix(SPECIESRAX.out.versions)

    // Now prepare for analysis with GeneRax
    ch_all_map_links = ch_core_og_maplinks.concat(ch_rem_og_maplinks)
    ch_all_gene_trees = ch_core_gene_trees.concat(ch_rem_gene_trees)
    ch_all_og_clean_msas = ch_core_og_clean_msas.concat(ch_rem_og_clean_msas)

    // Join these so that each gene family may be dealt with asynchronously as soon
    // as possible, and include with them the species tree.
    ch_generax_input = ch_all_map_links
        .join(ch_all_gene_trees)
        .join(ch_all_og_clean_msas)
        .combine(ch_speciesrax)

    // GENERAX_PER_FAMILY: full mode only
    if (params.generax_per_family) {
        GENERAX_PER_FAMILY(ch_generax_input)
        ch_versions = ch_versions.mix(GENERAX_PER_FAMILY.out.versions)
    }

    // GENERAX_PER_SPECIES: both modes
    GENERAX_PER_SPECIES(ch_generax_input)
    ch_versions = ch_versions.mix(GENERAX_PER_SPECIES.out.versions)

    ch_recon_perspp_gene_trees = GENERAX_PER_SPECIES.out.generax_per_spp_gfts.collect { it[1] }

    //
    // MODULE: PHYLO_PROFILES (batched)
    // Generate phylogenetic profiles from GeneRax reconciliation outputs
    // Batch inputs to avoid staging too many files at once
    //
    def batch_size = 1000

    // Combine all related data for each orthogroup into a tuple, then batch
    ch_phylo_profiles_input = GENERAX_PER_SPECIES.out.event_counts
        .join(GENERAX_PER_SPECIES.out.species_event_counts)
        .join(GENERAX_PER_SPECIES.out.transfer_event_counts)
        .join(GENERAX_PER_SPECIES.out.species_coverage)
        .map { meta, event_count, species_event_count, transfer_event_count, species_coverage ->
            [meta.og, event_count, species_event_count, transfer_event_count, species_coverage]
        }
        .toList()
        .flatMap { items ->
            items.collate(batch_size).withIndex().collect { batch, idx ->
                def ogs = batch.collect { it[0] }
                def event_counts = batch.collect { it[1] }
                def species_event_counts = batch.collect { it[2] }
                def transfer_event_counts = batch.collect { it[3] }
                def species_coverages = batch.collect { it[4] }
                [idx, event_counts, species_event_counts, transfer_event_counts, species_coverages, ogs]
            }
        }

    PHYLO_PROFILES(
        ch_phylo_profiles_input,
        ORTHOFINDER_MCL_ALL.out.inflation_dir
    )

    // Merge batched outputs
    MERGE_PHYLO_PROFILES(
        PHYLO_PROFILES.out.duplication_count.collect(),
        PHYLO_PROFILES.out.hgt_summed_count.collect(),
        PHYLO_PROFILES.out.loss_count.collect(),
        PHYLO_PROFILES.out.speciation_count.collect(),
        PHYLO_PROFILES.out.transfer_donor_count.collect(),
        PHYLO_PROFILES.out.transfer_recipient_count.collect()
    )
    ch_versions = ch_versions.mix(PHYLO_PROFILES.out.versions)
    ch_versions = ch_versions.mix(MERGE_PHYLO_PROFILES.out.versions)

    //
    // MODULE: PHYSICOCHEMICAL_PROPS
    // Calculate physicochemical properties for all gene families
    //
    if (params.zoogle) {
        PHYSICOCHEMICAL_PROPS(
            ch_all_og_clean_msas
        )
        ch_versions = ch_versions.mix(PHYSICOCHEMICAL_PROPS.out.versions)

        //
        // MODULE: BUILD_REFERENCE_CHRONOGRAM / TIME_CALIBRATE_SPECIES_TREE
        // Either auto-build a reference chronogram from TimeTree.org or use user-provided tree,
        // then time-calibrate the species tree for phylogenetic distance analysis
        //
        if (!params.reference_time_tree || params.reference_time_tree == 'none') {
            // Auto-build reference chronogram from TimeTree.org
            ch_species_names_file = species_name_list
                .collectFile(name: 'species_names.txt', newLine: true)
            BUILD_REFERENCE_CHRONOGRAM(ch_species_names_file, params.ncbi_email)
            ch_reference_tree = BUILD_REFERENCE_CHRONOGRAM.out.chronogram
            ch_versions = ch_versions.mix(BUILD_REFERENCE_CHRONOGRAM.out.versions)
        } else {
            // Use user-provided reference time tree
            ch_reference_tree = file(params.reference_time_tree)
        }

        TIME_CALIBRATE_SPECIES_TREE(
            ch_speciesrax,
            ch_reference_tree,
            params.time_calibration_method
        )
        ch_versions = ch_versions.mix(TIME_CALIBRATE_SPECIES_TREE.out.versions)

        //
        // MODULE: DATE_GENE_FAMILY_TREES
        // Time-calibrate gene family trees using reconciliation-filtered calibrations.
        // Only speciation nodes (S, SL) from GeneRax NHX trees are used;
        // duplication and transfer nodes are excluded.
        //
        // Build input: join NHX + Newick + MSA per OG
        ch_dating_input = GENERAX_PER_SPECIES.out.generax_per_spp_events       // [meta, nhx]
            .join(GENERAX_PER_SPECIES.out.generax_per_spp_gfts)             // [meta, nhx, newick]
            .join(ch_all_og_clean_msas)                                      // [meta, nhx, newick, msa]

        DATE_GENE_FAMILY_TREES(
            ch_dating_input,
            TIME_CALIBRATE_SPECIES_TREE.out.calibrated_tree,
            params.max_treepl_tips,
            params.age_bracket
        )
        ch_versions = ch_versions.mix(DATE_GENE_FAMILY_TREES.out.versions)

        // Feed dated trees into PHYLO_DIST (replaces raw reconciled trees)
        ch_phylo_dist_input = DATE_GENE_FAMILY_TREES.out.dated_gft          // [meta, dated_tree]
            .join(PHYSICOCHEMICAL_PROPS.out.summary_stats)                   // [meta, dated_tree, props]
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
                def refCount = proteinIds.count { it.startsWith("${params.ref_species}_") }

                // Count non-reference proteins
                def nonrefCount = proteinIds.size() - refCount

                // Count proteins per non-reference species
                // NOTE: Must match R script's species extraction logic (line 145 of protein_distance_calculation_functions.R)
                // Protein labels have format: Species_name_ProteinID
                // Species names are extracted by removing the last underscore-delimited segment
                def nonrefProteinsBySpecies = proteinIds
                    .findAll { !it.startsWith("${params.ref_species}_") }
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

        PHYLO_DIST(
            ch_phylo_dist_input,
            params.ref_species
        )
        ch_versions = ch_versions.mix(PHYLO_DIST.out.versions)
    }

    //
    // MODULE: ORTHOFINDER_PHYLOHOGS
    // Now using the reconciled gene family trees and rooted species tree,
    // parse orthogroups/gene families into hierarchical orthogroups (HOGs)
    // to identify orthologs and output orthogroup-level summary stats.
    //
    ORTHOFINDER_PHYLOHOGS(
        ch_speciesrax,
        ORTHOFINDER_MCL_ALL.out.inflation_dir,
        ORTHOFINDER_PREP_ALL.out.fastas,
        ORTHOFINDER_PREP_ALL.out.sppIDs,
        ORTHOFINDER_PREP_ALL.out.seqIDs,
        ch_recon_perspp_gene_trees,
        DIAMOND_BLASTP_ALL.out.txt.collect()
    )
}

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NOVELTREE()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
