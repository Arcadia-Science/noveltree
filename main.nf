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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT PARAMETER-SPECIFIED ALTERNATIVE MODULES (INCLUDES LOCAL AND NF-CORE-MODIFIED)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// TODO: Build into subworkflows
ch_aligner = params.aligner
ch_msa_trimmer = params.msa_trimmer
ch_tree_method = params.tree_method
if (ch_aligner == "witch") {
    include { WITCH as ALIGN_SEQS                   } from './modules/local/witch'
    include { WITCH as ALIGN_REMAINING_SEQS         } from './modules/local/witch'
} else {
    include { MAFFT as ALIGN_SEQS                   } from './modules/nf-core-modified/mafft'
    include { MAFFT as ALIGN_REMAINING_SEQS         } from './modules/nf-core-modified/mafft'
}
if (ch_msa_trimmer == "clipkit") {
    include { CLIPKIT as TRIM_MSAS                  } from './modules/local/clipkit'
    include { CLIPKIT as TRIM_REMAINING_MSAS        } from './modules/local/clipkit'
} else if (ch_msa_trimmer == 'cialign') {
    include { CIALIGN as TRIM_MSAS                  } from './modules/local/cialign'
    include { CIALIGN as TRIM_REMAINING_MSAS        } from './modules/local/cialign'
}
if (ch_tree_method == "iqtree") {
    include { IQTREE as INFER_TREES                 } from './modules/nf-core-modified/iqtree'
    include { IQTREE as INFER_REMAINING_TREES       } from './modules/nf-core-modified/iqtree'
} else {
    include { FASTTREE as INFER_TREES               } from './modules/local/fasttree'
    include { FASTTREE as INFER_REMAINING_TREES     } from './modules/local/fasttree'
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
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_ALL  } from './modules/local/orthofinder_prep'
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_TEST } from './modules/local/orthofinder_prep'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_TEST   } from './modules/local/orthofinder_mcl'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_ALL    } from './modules/local/orthofinder_mcl'
include { ANNOTATE_UNIPROT                          } from './modules/local/annotate_uniprot'
include { COGEQC                                    } from './modules/local/cogeqc'
include { SELECT_INFLATION                          } from './modules/local/select_inflation'
include { FILTER_ORTHOGROUPS                        } from './modules/local/filter_orthogroups'
include { ASTEROID                                  } from './modules/local/asteroid'
include { SPECIESRAX                                } from './modules/local/speciesrax'
include { GENERAX_PER_FAMILY                        } from './modules/local/generax_per_family'
include { GENERAX_PER_SPECIES                       } from './modules/local/generax_per_species'
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
include { DIAMOND_BLASTP as DIAMOND_BLASTP_ALL      } from './modules/nf-core-modified/diamond_blastp'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_TEST     } from './modules/nf-core-modified/diamond_blastp'
include { IQTREE_PMSF as IQTREE_PMSF_ALL         } from './modules/nf-core-modified/iqtree_pmsf'
include { IQTREE_PMSF as IQTREE_PMSF_REMAINING      } from './modules/nf-core-modified/iqtree_pmsf'

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
        ch_inflation = Channel.fromList(mcl_inflation)
        mcl_test_prots_list = ch_all_data.mcl_test_prots.collect { it[1] }

        ch_all_data.uniprot_prots.ifEmpty {
            error('Sample sheet must include samples from UniProt (uniprot column value is true) when doing MCL parameter best fit!')
        }

        //
        // MODULE: Annotate UniProt Proteins
        //
        ANNOTATE_UNIPROT(ch_all_data.uniprot_prots)
            .cogeqc_annotations
            .collect()
            .set { ch_annotations }
        ch_versions = ch_versions.mix(ANNOTATE_UNIPROT.out.versions)

        ORTHOFINDER_PREP_TEST(mcl_test_prots_list, "mcl_test_dataset")

        // Run for the test set (used to determine the best value of the MCL
        // inflation parameter)
        DIAMOND_BLASTP_TEST(
            ch_all_data.mcl_test_prots,
            ORTHOFINDER_PREP_TEST.out.fastas.flatten(),
            ORTHOFINDER_PREP_TEST.out.diamonds.flatten(),
            "txt",
            "true"
        )

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
    } else {
        ch_best_inflation = Channel.of(mcl_inflation.first())
    }

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
        params.min_num_grp_per_og,
        params.max_copy_num_spp_tree,
        params.max_copy_num_gene_trees
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
    // MODULE: ALIGN_SEQS
    // Infer multiple sequence alignments of orthogroups/gene
    // families using WITCH (default) or MAFFT
    //
    // For the extreme core set to be used in species tree inference
    if (ch_aligner == "witch") {
        ALIGN_SEQS(ch_spptree_fas)
        ch_versions = ch_versions.mix(ALIGN_SEQS.out.versions)
    } else {
        ALIGN_SEQS(ch_spptree_fas)
        ch_versions = ch_versions.mix(ALIGN_SEQS.out.versions)
    }

    // And for the remaining orthogroups:
    // Only start once species tree MSAs have finished (to give them priority)
    // We use the combination of collect().count() to hold off on running this
    // set of MSAs, while avoiding unnecessarily staging thousands of large files.
    if (ch_aligner == "witch") {
        ALIGN_REMAINING_SEQS(ch_genetree_fas)
    } else {
        ALIGN_REMAINING_SEQS(ch_genetree_fas)
    }

    //
    // MODULE: TRIM_MSAS
    // Trim gappy regions, poorly aligned, or and phylogenetically
    // uninformative/problematic sites from the MSAs using either
    // CIAlign or ClipKIT based on parameter specification.
    //
    if (ch_msa_trimmer == 'none') {
        if (ch_aligner == 'witch') {
            ch_core_og_maplinks = ALIGN_SEQS.out.map_link
            ch_rem_og_maplinks = ALIGN_REMAINING_SEQS.out.map_link
            ch_core_og_clean_msas = ALIGN_SEQS.out.cleaned_msas
            ch_rem_og_clean_msas = ALIGN_REMAINING_SEQS.out.cleaned_msas
        } else {
            ch_core_og_maplinks = ALIGN_SEQS.out.map_link
            ch_rem_og_maplinks = ALIGN_REMAINING_SEQS.out.map_link
            ch_core_og_clean_msas = ALIGN_SEQS.out.msas
            ch_rem_og_clean_msas = ALIGN_REMAINING_SEQS.out.msas
        }
    } else {
        if (ch_aligner == 'witch') {
            TRIM_MSAS(ALIGN_SEQS.out.cleaned_msas)
            TRIM_REMAINING_MSAS(ALIGN_REMAINING_SEQS.out.cleaned_msas)
            ch_core_og_maplinks = TRIM_MSAS.out.map_link
            ch_rem_og_maplinks = TRIM_REMAINING_MSAS.out.map_link
            ch_core_og_clean_msas = TRIM_MSAS.out.cleaned_msas
            ch_rem_og_clean_msas = TRIM_REMAINING_MSAS.out.cleaned_msas
            ch_versions = ch_versions.mix(TRIM_MSAS.out.versions)
        } else {
            TRIM_MSAS(ALIGN_SEQS.out.msas)
            TRIM_REMAINING_MSAS(ALIGN_REMAINING_SEQS.out.msas)
            ch_core_og_maplinks = TRIM_MSAS.out.map_link
            ch_rem_og_maplinks = TRIM_REMAINING_MSAS.out.map_link
            ch_core_og_clean_msas = TRIM_MSAS.out.cleaned_msas
            ch_rem_og_clean_msas = TRIM_REMAINING_MSAS.out.cleaned_msas
            ch_versions = ch_versions.mix(TRIM_MSAS.out.versions)
        }
    }
    // Create channels that are just lists of all the msas, and protein-species
    // map links that are provided in bulk to SpeciesRax
    core_og_maplink_list = ch_core_og_maplinks.collect { it[1] }
    core_og_clean_msa_list = ch_core_og_clean_msas.collect { it[1] }

    //
    // MODULE: INFER_TREES
    // Infer gene-family trees from the trimmed MSAs using either
    // VeryFastTree or IQ-TREE.
    //
    INFER_TREES(ch_core_og_clean_msas, params.tree_model)
    INFER_REMAINING_TREES(ch_rem_og_clean_msas, params.tree_model)
    ch_versions = ch_versions.mix(INFER_TREES.out.versions)

    // Run IQ-TREE PMSF if model is specified, and subsequently collect final
    // phylogenies into a channel for downstram use
    if (params.tree_model_pmsf != 'none') {
        //
        // MODULE: IQTREE_PMSF
        // Infer gene-family trees from the trimmed MSAs and guide trees from the
        // previous tree inference module
        //
        // Be sure that both the MSAs and guide trees are sorted into the same
        // order as before to prevent any hiccups - do so by temporarily
        // joining the two channels.
        ch_pmsf_input = ch_core_og_clean_msas.join(INFER_TREES.out.phylogeny)
        ch_pmsf_input_remaining = ch_rem_og_clean_msas.join(INFER_REMAINING_TREES.out.phylogeny)
        // Now run
        IQTREE_PMSF(ch_pmsf_input, params.tree_model_pmsf)
        IQTREE_PMSF_REMAINING(ch_pmsf_input_remaining, params.tree_model_pmsf)
        ch_versions = ch_versions.mix(IQTREE_PMSF.out.versions)

        ch_core_gene_trees = IQTREE_PMSF.out.phylogeny
        ch_rem_gene_trees = IQTREE_PMSF_REMAINING.out.phylogeny
        // And create a channel/list (no tuple) of just the core trees used by Asteroid
        core_gene_tree_list = ch_core_gene_trees.collect { it[1] }
    } else {
        ch_core_gene_trees = INFER_TREES.out.phylogeny
        ch_rem_gene_trees = INFER_REMAINING_TREES.out.phylogeny
        core_gene_tree_list = ch_core_gene_trees.collect { it[1] }
    }

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
    ch_all_map_links = ch_core_og_maplinks
        .concat(ch_rem_og_maplinks)
    ch_all_gene_trees = ch_core_gene_trees
        .concat(ch_rem_gene_trees)
    ch_all_og_clean_msas = ch_core_og_clean_msas
        .concat(ch_rem_og_clean_msas)

    // Join these so that each gene family may be dealt with asynchronously as soon
    // as possible, and include with them the species tree.
    ch_generax_input = ch_all_map_links
        .join(ch_all_gene_trees)
        .join(ch_all_og_clean_msas)
        .combine(ch_speciesrax)

    GENERAX_PER_FAMILY(
        ch_generax_input
    )
        .generax_per_fam_gfts
        .collect { it[1] }
        .set { ch_recon_perfam_gene_trees }

    GENERAX_PER_SPECIES(
        ch_generax_input
    )
        .generax_per_spp_gfts
        .collect { it[1] }
        .set { ch_recon_perspp_gene_trees }

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
