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
if ((params.min_prop_spp_for_spptree as double) <= 0.0 ||
    (params.min_prop_spp_for_spptree as double) > 1.0) {
    exit 1, '--min_prop_spp_for_spptree must be in (0, 1]'
}
if ((params.max_copy_num_spp_tree as double) < 1.0) {
    exit 1, '--max_copy_num_spp_tree must be at least 1'
}
if ((params.speciesrax_min_species_occupancy as double) <= 0.0 ||
    (params.speciesrax_min_species_occupancy as double) > 1.0) {
    exit 1, '--speciesrax_min_species_occupancy must be in (0, 1]'
}
if ((params.speciesrax_max_mean_copies as double) < 1.0 ||
    (params.speciesrax_max_copies_per_species as int) < 1 ||
    (params.speciesrax_max_total_leaves_factor as double) < 1.0) {
    exit 1, 'SpeciesRax copy-number and total-leaf limits must be at least 1'
}
if ((params.speciesrax_bundle_size as int) < 1) {
    exit 1, '--speciesrax_bundle_size must be at least 1'
}
if ((params.align_tier1_max as int) < 1 ||
    (params.align_tier2_max as int) < (params.align_tier1_max as int)) {
    exit 1, 'Adaptive alignment thresholds must satisfy 1 <= --align_tier1_max <= --align_tier2_max'
}
// A trusted root is required before SpeciesRax. Zoogle also reuses this exact
// chronogram for dating, so root inference and calibration cannot disagree.
def needs_reference_chronogram = params.zoogle || params.outgroups == 'none'
if (needs_reference_chronogram && (!params.reference_time_tree || params.reference_time_tree == 'none')) {
    if (!params.ncbi_email || params.ncbi_email == 'none') {
        exit 1, 'Species-tree rooting requires --outgroups, --reference_time_tree, or --ncbi_email to auto-build a rooted TimeTree chronogram. SpeciesRax no longer guesses the root with --si-strategy REROOT.'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS AND MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK                                    } from './subworkflows/local/input_check'
include { PREPARE_INPUTS                                 } from './subworkflows/local/prepare_inputs'
include { INFER_ORTHOGROUPS                              } from './subworkflows/local/infer_orthogroups'
include { INFER_GENE_TREES as SPECIESTREE_GENE_FAMILIES  } from './subworkflows/local/infer_gene_trees'
include { INFER_GENE_TREES as REMAINING_GENE_FAMILIES    } from './subworkflows/local/infer_gene_trees'
include { RECONCILE_TREES                                } from './subworkflows/local/reconcile_trees'
include { RECONCILIATION_SUMMARIES                       } from './subworkflows/local/reconciliation_summaries'
include { BUILD_REFERENCE_CHRONOGRAM                     } from './modules/local/build_reference_chronogram'

if (params.zoogle) {
    include { ZOOGLE } from './subworkflows/local/zoogle'
}

if (params.busco) {
    include { BUSCO as BUSCO_SHALLOW } from './modules/nf-core-modified/busco'
    include { BUSCO as BUSCO_BROAD   } from './modules/nf-core-modified/busco'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NOVELTREE {

    // 1. Validate samplesheet
    ch_all_data = INPUT_CHECK(ch_input)

    // Normalize ref_species to hyphens (users may pass underscores or spaces)
    // Set to 'none' for centroid-only analysis (no reference species)
    def ref_species = (params.ref_species && params.ref_species != 'none')
        ? params.ref_species.replace('_', '-').replace(' ', '-')
        : 'none'

    // 2. Download, preprocess, rename
    PREPARE_INPUTS(ch_all_data.remote_prots, ch_all_data.local_prots)
    ch_renamed_prots = PREPARE_INPUTS.out.renamed_prots
    ch_protein_maps = PREPARE_INPUTS.out.protein_maps

    // These value channels are consumed by multiple downstream subworkflows
    species_name_list   = ch_renamed_prots.collect { it[0].id }
    complete_prots_list = ch_renamed_prots.collect { it[1] }

    // Prepare the trusted chronogram once, before reconciliation. When no
    // explicit outgroup is supplied, its encoded root split roots the
    // MiniNJ topology passed to the final SpeciesRax fit. Zoogle later uses the
    // same file for calibration ages.
    if (needs_reference_chronogram) {
        if (!params.reference_time_tree || params.reference_time_tree == 'none') {
            ch_species_names_file = species_name_list
                .collectFile(name: 'species_names.txt', newLine: true)
            BUILD_REFERENCE_CHRONOGRAM(ch_species_names_file, params.ncbi_email)
            ch_reference_tree = BUILD_REFERENCE_CHRONOGRAM.out.chronogram
        } else {
            ch_reference_tree = Channel.value(file(params.reference_time_tree, checkIfExists: true))
        }
    } else {
        // An explicit outgroup roots SpeciesRax without a reference file.
        // Pass an empty optional path rather than staging a fake sentinel.
        ch_reference_tree = Channel.value([])
    }

    // Warn if ref_species is set but not found in the dataset
    if (params.zoogle && ref_species != 'none') {
        species_name_list.collect().subscribe { names ->
            if (!names.contains(ref_species)) {
                log.warn "WARNING: ref_species '${ref_species}' not found in input dataset. " +
                         "Reference-based analyses will be skipped. Only centroid-based outputs will be produced."
            }
        }
    }

    // 3. BUSCO (optional, independent — no downstream consumers)
    if (params.busco) {
        BUSCO_SHALLOW(
            ch_renamed_prots.filter{ it[0].busco_shallow != "NA" },
            "shallow",
            [],
            []
        )
        BUSCO_BROAD(
            ch_renamed_prots.filter{ it[0].busco_broad != "NA" },
            "broad",
            [],
            []
        )
    }

    // 4. Orthogroup inference (--test_run filters to all-species OGs only)
    if (params.test_run) {
        log.info "── TEST MODE: only gene families containing ALL species will be processed ──"
    }
    INFER_ORTHOGROUPS(
        ch_renamed_prots,
        complete_prots_list,
        mcl_inflation,
        INPUT_CHECK.out.complete_samplesheet,
        ch_protein_maps.collect { it[1] }
    )

    // 5. Gene tree inference (alignment → trimming → phylogeny)
    SPECIESTREE_GENE_FAMILIES(INFER_ORTHOGROUPS.out.spptree_fas)
    REMAINING_GENE_FAMILIES(INFER_ORTHOGROUPS.out.genetree_fas)

    // 6. Species tree + gene-tree/species-tree reconciliation
    RECONCILE_TREES(
        species_name_list,
        SPECIESTREE_GENE_FAMILIES.out.phylogeny,
        SPECIESTREE_GENE_FAMILIES.out.map_link,
        SPECIESTREE_GENE_FAMILIES.out.cleaned_msas,
        REMAINING_GENE_FAMILIES.out.phylogeny,
        REMAINING_GENE_FAMILIES.out.map_link,
        REMAINING_GENE_FAMILIES.out.cleaned_msas,
        ch_reference_tree,
        INFER_ORTHOGROUPS.out.speciesrax_selection,
        INFER_ORTHOGROUPS.out.speciesrax_coverage
    )

    // 7. Reconciliation summaries (phylo profiles, HOG parsing)
    RECONCILIATION_SUMMARIES(
        RECONCILE_TREES.out.species_event_counts,
        RECONCILE_TREES.out.species_coverage,
        RECONCILE_TREES.out.generax_nhx,
        RECONCILE_TREES.out.labeled_species_tree
    )

    // 8. Zoogle: physicochemical properties + phylogenetic distance analysis (conditional)
    if (params.zoogle) {
        ZOOGLE(
            INFER_ORTHOGROUPS.out.spptree_fas,
            INFER_ORTHOGROUPS.out.genetree_fas,
            RECONCILE_TREES.out.all_clean_msas,
            RECONCILE_TREES.out.generax_per_spp_gfts,
            RECONCILE_TREES.out.speciesrax_tree,
            species_name_list,
            ref_species,
            RECONCILIATION_SUMMARIES.out.orthologs,
            RECONCILIATION_SUMMARIES.out.paralogs,
            ch_reference_tree
        )
    }
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
