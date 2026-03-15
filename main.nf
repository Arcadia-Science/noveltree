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
    def ref_species = params.ref_species.replace('_', '-').replace(' ', '-')

    // 2. Download, preprocess, rename
    PREPARE_INPUTS(ch_all_data.remote_prots, ch_all_data.local_prots)
    ch_renamed_prots = PREPARE_INPUTS.out.renamed_prots

    // These value channels are consumed by multiple downstream subworkflows
    species_name_list   = ch_renamed_prots.collect { it[0].id }
    complete_prots_list = ch_renamed_prots.collect { it[1] }

    // 3. BUSCO (optional, independent — no downstream consumers)
    if (params.busco) {
        BUSCO_SHALLOW(
            ch_renamed_prots.filter{ it[0].shallow_db != "NA" },
            "shallow",
            [],
            []
        )
        BUSCO_BROAD(
            ch_renamed_prots.filter{ it[0].broad_db != "NA" },
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
        INPUT_CHECK.out.complete_samplesheet
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
        REMAINING_GENE_FAMILIES.out.cleaned_msas
    )

    // 7. Reconciliation summaries (phylo profiles, HGT matrix, HOG parsing)
    RECONCILIATION_SUMMARIES(
        RECONCILE_TREES.out.event_counts,
        RECONCILE_TREES.out.species_event_counts,
        RECONCILE_TREES.out.transfer_event_counts,
        RECONCILE_TREES.out.species_coverage,
        INFER_ORTHOGROUPS.out.inflation_dir,
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
            ref_species
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
