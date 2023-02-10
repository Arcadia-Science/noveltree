#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Arcadia-Science/phylorthology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Arcadia-Science/phylorthology
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
    ch_mcl_inflation = Channel.of(params.mcl_inflation.split(","))
} else {
    exit 1, 'MCL Inflation parameter(s) not specified!'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW
//
include { INPUT_CHECK                               } from './subworkflows/local/input_check'

//
// MODULE
//
// Modules being run twice (for MCL testing and full analysis)
// needs to be included twice under different names.
include { ORTHOFINDER_PREP                          } from './modules/local/orthofinder_prep'
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_TEST } from './modules/local/orthofinder_prep'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_TEST   } from './modules/local/orthofinder_mcl'
include { ORTHOFINDER_MCL                           } from './modules/local/orthofinder_mcl'
include { ANNOTATE_UNIPROT                          } from './modules/local/annotate_uniprot'
include { COGEQC                                    } from './modules/local/cogeqc'
include { SELECT_INFLATION                          } from './modules/local/select_inflation'
include { FILTER_ORTHOGROUPS                        } from './modules/local/filter_orthogroups'
include { CLIPKIT                                   } from './modules/local/clipkit'
include { CLIPKIT as CLIPKIT_REMAINING              } from './modules/local/clipkit'
include { SPECIES_TREE_PREP                         } from './modules/local/species_tree_prep'
include { SPECIES_TREE_PREP as GENE_TREE_PREP       } from './modules/local/species_tree_prep'
include { ASTEROID                                  } from './modules/local/asteroid'
include { SPECIESRAX                                } from './modules/local/speciesrax'
include { GENERAX                                   } from './modules/local/generax'
include { ORTHOFINDER_PHYLOHOGS                     } from './modules/local/orthofinder_phylohogs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE-MODIFIED MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE
//
// Modules being run twice (for MCL testing and full analysis)
// needs to be included twice under different names.
include { BUSCO as BUSCO_SHALLOW                    } from './modules/nf-core-modified/busco'
include { BUSCO as BUSCO_BROAD                      } from './modules/nf-core-modified/busco'
include { DIAMOND_BLASTP                            } from './modules/nf-core-modified/diamond_blastp'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_TEST     } from './modules/nf-core-modified/diamond_blastp'
include { IQTREE                                    } from './modules/nf-core-modified/iqtree'
include { IQTREE as IQTREE_REMAINING                } from './modules/nf-core-modified/iqtree'
include { MAFFT                                     } from './modules/nf-core-modified/mafft'
include { MAFFT as MAFFT_REMAINING                  } from './modules/nf-core-modified/mafft'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main Arcadia-Science/phylorthology analysis pipeline
//
workflow PHYLORTHOLOGY {
    ch_inflation = ch_mcl_inflation.toList().flatten()
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_all_data = INPUT_CHECK(ch_input)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    complete_prots_list = ch_all_data.complete_prots.collect { it[1] }
    mcl_test_prots_list = ch_all_data.mcl_test_prots.collect { it[1] }
    uniprot_prots_list = ch_all_data.uniprot_prots.collect { it[1] }

    //
    // MODULE: Run BUSCO
    // Split up into shallow and broad scale runs, since downstream modules
    // do not use these outputs, so multiple busco runs may be conducted
    // simultaneously
    //
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

    //
    // MODULE: Annotate UniProt Proteins
    //
    ANNOTATE_UNIPROT(ch_all_data.uniprot_prots, params.download_annots)
        .cogeqc_annotations
        .collect()
        .set { ch_annotations }
    ch_versions = ch_versions.mix(ANNOTATE_UNIPROT.out.versions)

    //
    // MODULE: Prepare directory structure and fasta files according to
    //         OrthoFinder's preferred format for downstream MCL clustering
    //
    ORTHOFINDER_PREP(complete_prots_list, "complete_dataset")
    ORTHOFINDER_PREP_TEST(mcl_test_prots_list, "mcl_test_dataset")
    ch_versions = ch_versions.mix(ORTHOFINDER_PREP.out.versions)

    //
    // MODULE: All-v-All diamond/blastp
    //
    // Run for the test set (used to determine the best value of the MCL
    // inflation parameter)
    DIAMOND_BLASTP_TEST(
        ch_all_data.mcl_test_prots,
        ORTHOFINDER_PREP_TEST.out.fastas.flatten(),
        ORTHOFINDER_PREP_TEST.out.diamonds.flatten(),
        "txt",
        "true"
    )

    // And for the full dataset, to be clustered into orthogroups using
    // the best inflation parameter.
    DIAMOND_BLASTP(
        ch_all_data.complete_prots,
        ORTHOFINDER_PREP.out.fastas.flatten(),
        ORTHOFINDER_PREP.out.diamonds.flatten(),
        "txt",
        "false"
    )
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

    //
    // MODULE: Run Orthofinder's implementation of MCL (with similarity score
    //         correction).
    //

    // TODO: Fix the output_dir determination logic
    // First determine the optimal MCL inflation parameter, and then
    // subsequently use this for full orthogroup inference.
    ORTHOFINDER_MCL_TEST(
        ch_inflation,
        DIAMOND_BLASTP_TEST.out.txt.collect(),
        ORTHOFINDER_PREP_TEST.out.fastas,
        ORTHOFINDER_PREP_TEST.out.diamonds,
        ORTHOFINDER_PREP_TEST.out.sppIDs,
        ORTHOFINDER_PREP_TEST.out.seqIDs,
        "mcl_test_dataset"
    )

    //
    // MODULE: COGEQC
    // Run an R-script that applies cogqc to assess orthogroup inference
    // accuracy/performance.
    //
    COGEQC(
        ORTHOFINDER_MCL_TEST.out.inflation_dir,
        params.min_num_spp_per_og,
        ch_annotations
    )
    ch_cogeqc_summary = COGEQC.out.cogeqc_summary.collect()
    ch_versions = ch_versions.mix(COGEQC.out.versions)

    // Now, from these orthogroup summaries, select the best inflation parameter
    SELECT_INFLATION(ch_cogeqc_summary, params.min_num_spp_per_og)
        .best_inflation.text.trim()
        .set { ch_best_inflation }
    ch_versions = ch_versions.mix(SELECT_INFLATION.out.versions)

    // Using this best-performing inflation parameter, infer orthogroups for
    // all samples.
    ORTHOFINDER_MCL(
        ch_best_inflation,
        DIAMOND_BLASTP.out.txt.collect(),
        ORTHOFINDER_PREP.out.fastas,
        ORTHOFINDER_PREP.out.diamonds,
        ORTHOFINDER_PREP.out.sppIDs,
        ORTHOFINDER_PREP.out.seqIDs,
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
        ORTHOFINDER_MCL.out.inflation_dir,
        params.min_num_spp_per_og,
        params.min_num_grp_per_og,
        params.max_copy_num_spp_tree,
        params.max_copy_num_gene_trees
    )

    //
    // MODULE: MAFFT
    // Infer multiple sequence alignments of orthogroups/gene
    // families using MAFFT
    //
    // For the extreme core set to be used in species tree inference
    MAFFT(
        FILTER_ORTHOGROUPS.out.spptree_fas.flatten(),
        []
    )
        .msas
        .set{ ch_core_og_msas }

    // And for the remaining orthogroups:
    // Only start once species tree MSAs have finished (to give them priority)
    // We use the combination of collect().count() to hold off on running this 
    // set of MSAs, while avoiding unnecessarily staging thousands of large files.
    MAFFT_REMAINING(
        FILTER_ORTHOGROUPS.out.genetree_fas.flatten(),
        ch_core_og_msas.collect().count()
    )
        .msas
        .set { ch_rem_og_msas }
    ch_versions = ch_versions.mix(MAFFT.out.versions)

    //
    //MODULE: CLIPKIT
    // Trim gappy and phylogenetically uninformative sites from the MSAs
    //
    ch_core_trimmed_msas = CLIPKIT(ch_core_og_msas).trimmed_msas

    ch_rem_trimmed_msas = CLIPKIT_REMAINING(ch_rem_og_msas).trimmed_msas
    ch_versions = ch_versions.mix(CLIPKIT.out.versions)

    //
    // MODULE: IQTREE
    // Infer gene-family trees from the trimmed MSAs
    //
    IQTREE(
        ch_core_trimmed_msas,
        params.tree_model,
        params.tree_model_pmsf
    )
        .phylogeny
        .collect()
        .set { ch_core_gene_trees }

    IQTREE_REMAINING(
        ch_rem_trimmed_msas,
        params.tree_model,
        params.tree_model_pmsf
    )
        .phylogeny
        .collect()
        .set { ch_rem_gene_trees }

    ch_versions = ch_versions.mix(IQTREE.out.versions)


    // Now, go ahead and prepare input files for initial unrooted species
    // tree inference with Asteroid, rooted species-tree inference with
    // SpeciesRax, and gene-tree species-tree reconciliation and estimation
    // of gene family duplication transfer and loss with GeneRax.

    // Do this for both the core and non-core gene families.
    // All outputs are needed for species tree inference, but not for the
    // remainder.
    SPECIES_TREE_PREP(
        ch_core_gene_trees,
        ch_core_trimmed_msas.collect()
    )
        .set { ch_core_spptree_prep }

    ch_core_treefile = ch_core_spptree_prep.treefile
    ch_core_families = ch_core_spptree_prep.families
    ch_core_generax_map = ch_core_spptree_prep.generax_map
    ch_asteroid_map = ch_core_spptree_prep.asteroid_map

    GENE_TREE_PREP(
        ch_rem_gene_trees,
        ch_rem_trimmed_msas.collect()
    )
        .set { ch_rem_genetree_prep }

    ch_rem_treefile = ch_rem_genetree_prep.treefile
    ch_rem_families = ch_rem_genetree_prep.families
    ch_rem_generax_map = ch_rem_genetree_prep.generax_map

    // The following two steps will just be done for the core set of
    // orthogroups that will be used to infer the species tree
    //
    // MODULE: ASTEROID
    // Alrighty, now let's infer an intial, unrooted species tree using Asteroid
    //
    ASTEROID(
        ch_core_treefile,
        ch_asteroid_map
    )
        .spp_tree
        .set { ch_asteroid }
    ch_versions = ch_versions.mix(ASTEROID.out.versions)

    //
    // MODULE: SPECIESRAX
    // Now infer the rooted species tree with SpeciesRax,
    // reconcile gene family trees, and infer per-family
    // rates of gene-family duplication, transfer, and loss
    //
    SPECIESRAX(
        ch_asteroid,
        ch_core_generax_map,
        ch_core_gene_trees.collect(),
        ch_core_trimmed_msas.collect(),
        ch_core_families
    )
        .speciesrax_tree
        .set { ch_speciesrax }
    ch_versions = ch_versions.mix(SPECIESRAX.out.versions)

    // Run again, but this time only using the GeneRax component,
    // reconciling gene family trees with the rooted species tree
    // inferred from SpeciesRax for all remaining gene families
    GENERAX(
        ch_speciesrax,
        ch_rem_generax_map,
        ch_rem_gene_trees.collect(),
        ch_rem_trimmed_msas.collect(),
        ch_rem_families
    )

    //
    // MODULE: ORTHOFINDER_PHYLOHOGS
    // Now using the reconciled gene family trees and rooted species tree,
    // parse orthogroups/gene families into hierarchical orthogroups (HOGs)
    // to identify orthologs and output orthogroup-level summary stats.
    //
    ORTHOFINDER_PHYLOHOGS(
        ch_speciesrax,
        ORTHOFINDER_MCL.out.inflation_dir,
        ORTHOFINDER_PREP.out.fastas,
        ORTHOFINDER_PREP.out.sppIDs,
        ORTHOFINDER_PREP.out.seqIDs,
        SPECIESRAX.out.speciesrax_gfts,
        GENERAX.out.generax_gfts,
        DIAMOND_BLASTP.out.txt.collect()
    )
}

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    PHYLORTHOLOGY()
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
