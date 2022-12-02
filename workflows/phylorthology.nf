/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowPhylorthology.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta_dir) { ch_fa_dir = params.fasta_dir } else { exit 1, 'Fasta directory not specified!' }
if (params.test_fasta_dir) { ch_test_fa_dir = params.test_fasta_dir } else { exit 1, 'Test fasta directory not specified!' }
if (params.mcl_inflation) { ch_mcl_inflation = Channel.of(params.mcl_inflation) } else { exit 1, 'MCL Inflation parameter(s) not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW
//
include { INPUT_CHECK                       } from '../subworkflows/local/input_check'

//
// MODULE
//
// Modules breing run twice (for MCL testing and fulla analysis)
// needs to be included twice under different names.
include { BUSCO as BUSCO_SHALLOW                    } from '../modules/local/busco'
include { BUSCO as BUSCO_BROAD                      } from '../modules/local/busco'
include { ORTHOFINDER_PREP                          } from '../modules/local/orthofinder_prep'
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_TEST } from '../modules/local/orthofinder_prep'
include { DIAMOND_BLASTP                            } from '../modules/local/diamond_blastp'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_TEST     } from '../modules/local/diamond_blastp'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_TEST   } from '../modules/local/orthofinder_mcl'
include { ORTHOFINDER_MCL                           } from '../modules/local/orthofinder_mcl'
include { ANNOTATE_UNIPROT                          } from '../modules/local/annotate_uniprot'
include { COGEQC                                    } from '../modules/local/cogeqc'
include { SELECT_INFLATION                          } from '../modules/local/select_inflation'
include { FILTER_ORTHOGROUPS                        } from '../modules/local/filter_orthogroups'
include { CLIPKIT                                   } from '../modules/local/clipkit'
include { CLIPKIT as CLIPKIT_REMAINING              } from '../modules/local/clipkit'
include { IQTREE                                    } from '../modules/local/iqtree'
include { IQTREE as IQTREE_REMAINING                } from '../modules/local/iqtree'
include { SPECIES_TREE_PREP                         } from '../modules/local/species_tree_prep'
include { SPECIES_TREE_PREP as GENE_TREE_PREP       } from '../modules/local/species_tree_prep'
include { ASTEROID                                  } from '../modules/local/asteroid'
include { SPECIESRAX                                } from '../modules/local/speciesrax'
include { GENERAX                                   } from '../modules/local/generax'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MAFFT                                 } from '../modules/nf-core/mafft/main' // Slightly modified (specific parameters)
include { MAFFT as MAFFT_REMAINING              } from '../modules/nf-core/mafft/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to get list of [ meta, [ file ] ]
// This function is applied to a csv produced by the "og-tax-summary.R"
// script called as part of the "filter_orthogroups.nf" local module.
def create_og_channel(LinkedHashMap row) {
    // create meta map
    def meta  = [:]
        meta.og   = row.orthogroup
        meta.nspp = row.num_spp
        meta.total_copy_num = row.total_copy_num
        meta.copy_num = row.mean_copy_num
        meta.num_grps = row.num_tax_grps
    // add path(s) of the OG file to the meta map
    def og_meta = []
        og_meta = [ meta, [ file(row.file) ] ] 
    return og_meta
}

workflow PHYLORTHOLOGY {
    
    ch_inflation = ch_mcl_inflation.toList().flatten()
    ch_versions = Channel.empty()
    
    //
    // MODULE: Prepare directory structure and fasta files according to 
    //         OrthoFinder's preferred format for downstream MCL clustering
    //
    ORTHOFINDER_PREP (
        ch_fa_dir,
        "false"
    )
    
    ORTHOFINDER_PREP_TEST (
        ch_test_fa_dir,
        "true"
    )
    
    // Fasta files should be redirected into a channel set of filepaths emitted
    // separately, whereas the diamond databases for each species can be put
    // into a directory as they are now (as a comma seperated list emitted 
    // together).
    ch_fa = ORTHOFINDER_PREP.out.fa.flatten()
    ch_dmd = ORTHOFINDER_PREP.out.dmd.flatten()
    ch_test_fa = ORTHOFINDER_PREP_TEST.out.fa.flatten()
    ch_test_dmd = ORTHOFINDER_PREP_TEST.out.dmd.flatten()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    // First the full set
    ch_prots_all = INPUT_CHECK (
        ch_input
    )
    .prots
    .map {
        meta, prots ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, prots ]
    }
    .groupTuple(by: [0])

    // Pull out the test and full sets
    ch_prots_all
    .branch {
        meta, prots ->
            proteomes  : meta.mcl_test == 'false'
                return [ meta, prots.flatten() ]
    }
    .set { ch_prots_full }
    
    ch_prots_all
    .branch {
        meta, prots ->
            proteomes  : meta.mcl_test == 'true'
                return [ meta, prots.flatten() ]
    }
    .set { ch_prots_test }
    
    ch_prots_full
    .branch {
        meta, prots ->
            proteomes  : prots
                return [ meta, prots.flatten() ]
    }
    .set { ch_busco_full }

    ch_prots_test
    .branch {
        meta, prots ->
            proteomes  : prots
                return [ meta, prots.flatten() ]
    }
    .set { ch_busco_test }
    
    // Create an orthofinder channel with paths to the new fasta/diamond DBs
    ch_orthofinder_full = ch_busco_full.merge(ch_fa, ch_dmd)
    ch_orthofinder_test = ch_busco_test.merge(ch_test_fa, ch_test_dmd)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Now, there will be a couple of modules below that we reapply, both to 
    // the full dataset, and to the MCL inflation parameter test set. 
    // These repeat modules include:
    // uniprot annotation (ch_annotations)
    // diamond blastp (ch_blastp & ch_similarities)
    
    //
    // MODULE: Annotate UniProt Proteins
    //
    ch_annotations = ANNOTATE_UNIPROT (
            ch_busco_full
        )
        .annotations
        .collect()
    ch_versions = ch_versions.mix(ANNOTATE_UNIPROT.out.versions)
    
    //
    // MODULE: Run BUSCO
    // Split up into shallow and broad scale runs, since downstream modules
    // do not use these outputs, so multiple busco runs may be conducted
    // simultaneously
    //
    // Shallow taxonomic scale:
    BUSCO_SHALLOW (
        ch_busco_full,
        "shallow",
        [],
        []
        )
    
    // Broad taxonomic scale (Eukaryotes)
    BUSCO_BROAD (
        ch_busco_full,
        "broad",
        [],
        []
        )
        
    //
    // MODULE: All-v-All diamond/blastp
    //
    // Run for the test set (used to determine the best value of the MCL
    // inflation parameter)
    ch_blastp_test = DIAMOND_BLASTP_TEST (
        ch_orthofinder_test,
        ch_test_dmd,
        "txt",
        "true",
        []
    )
    .txt
    
    // And for the full dataset, to be clustered into orthogroups using 
    // the best inflation parameter. 
    ch_blastp = DIAMOND_BLASTP (
        ch_orthofinder_full,
        ch_dmd,
        "txt",
        "false",
        []
    )
    .txt

    //
    // MODULE: Run Orthofinder's implementation of MCL (with similarity score
    //         correction).
    //
    // Collect all pairwise similarity scores into a single channel and pass to
    // the orthofinder MCL analysis so that it doesn't start until the full 
    // set of all-v-all comparisons have completed.
    ch_similarities_test = ch_blastp_test.mix(ch_blastp_test).collect()
    ch_similarities = ch_blastp.mix(ch_blastp).collect()

    // First determine the optimal MCL inflation parameter, and then 
    // subsequently use this for full orthogroup inference. 
    ch_mcl = ORTHOFINDER_MCL_TEST (
        ch_inflation,
        ch_similarities_test,
        "true"
    )
    .og_fpath
    
    //
    // MODULE: COGEQC
    // Run an R-script that applies cogqc to assess orthogroup inference
    // accuracy/performance.
    //
    ch_cogeqc = COGEQC (
        ch_mcl,
        ch_annotations
    )

    ch_summs = COGEQC.out.og_summary.collect()
    
    // Now, from these orthogroup summaries, select the best inflation parameter
    SELECT_INFLATION (
        ch_summs
    )
    .best_inflation
    .map{ file -> file.text.trim() } 
    .set { ch_best_inflation } 

    // Using this best-performing inflation parameter, infer orthogroups for
    // all samples. 
    ch_orthogroups = ORTHOFINDER_MCL (
        ch_best_inflation,
        ch_similarities,
        "false",
        )
        .og_fpath
    
    // 
    // MODULE: FILTER_ORTHOGROUPS
    // Subset orthogroups based on their copy number and distribution 
    // across species and taxonomic group. 
    // The conservative subset will be used for species tree inference, 
    // and the remainder will be used to infer gene family trees only. 
    ch_filtered_ogs = FILTER_ORTHOGROUPS (
        ch_orthogroups,
        "4",
        "4",
        "1",
        "2"
        )
    // Subset, pulling out two orthogroup sets:
    // one for species tree inference (core) and a remaining core set 
    // that we will infer gene family trees for (remaining (rem)). 
    // All 'core' gene family trees will be reconciled with the species tree, 
    // and duplication/tranfer/loss rates will be estimated for these, 
    // but not all orthogroups will have MSAs/gene family trees estimated 
    // (because they are either very taxon specific, or incredibly large, e.g. 
    // a mean per-species gene-copy number > 10).
    ch_filtered_ogs
    .spptree_core_ogs
    .splitCsv ( header:true, sep:',' )
    .map { create_og_channel(it) }
    .set { ch_core_ogs }
    
    ch_filtered_ogs
    .genetree_core_ogs
    .splitCsv ( header:true, sep:',' )
    .map { create_og_channel(it) }
    .set { ch_rem_ogs }
 
    //
    // MODULE: MAFFT
    // Infer multiple sequence alignments of orthogroups/gene 
    // families using MAFFT 
    //
    // For the extreme core set to be used in species tree inference
    ch_core_og_msas = MAFFT (
        ch_core_ogs
    )
    .msas
    
    // And for the remaining orthogroups
    ch_rem_og_msas = MAFFT_REMAINING (
        ch_rem_ogs
    )
    .msas
    
    //
    //MODULE: CLIPKIT 
    // Trim gappy and phylogenetically uninformative sites from the MSAs
    //
    ch_core_trimmed_msas = CLIPKIT (
        ch_core_og_msas
    )
    .trimmed_msas
    
    ch_rem_trimmed_msas = CLIPKIT_REMAINING (
        ch_rem_og_msas
    )
    .trimmed_msas
    
    //
    // MODULE: IQTREE
    // Infer gene-family trees from the trimmed MSAs
    //
    ch_core_gene_trees = IQTREE (
        ch_core_trimmed_msas,
        []
    )
    .phylogeny
    
    ch_rem_gene_trees = IQTREE_REMAINING (
        ch_rem_trimmed_msas,
        []
    )
    .phylogeny
    
    // Collect these gene family trees and alignments;
    // they will be used for unrooted species tree inference 
    // with Asteroid and downstream analysis with GeneRax and 
    // SpeciesRax
    
    // Do this for both the core and remaining orthogroups
    // First trees....
    ch_core_gene_trees
    .branch {
        meta, phylogeny ->
            trees  : phylogeny
                return phylogeny
    }
    .collect()
    .set { ch_all_core_trees } 
    
    ch_rem_gene_trees
    .branch {
        meta, phylogeny ->
            trees  : phylogeny
                return phylogeny
    }
    .collect()
    .set { ch_all_rem_trees } 
    
    // Then the alignments.
    ch_core_trimmed_msas
    .branch {
        meta, trimmed_msas ->
            msas  : trimmed_msas
                return trimmed_msas
    }
    .collect()
    .set { ch_all_core_msas } 
    
    ch_rem_trimmed_msas
    .branch {
        meta, trimmed_msas ->
            msas  : trimmed_msas
                return trimmed_msas
    }
    .collect()
    .set { ch_all_rem_msas }
    
    // Now, go ahead and prepare input files for initial unrooted species 
    // tree inference with Asteroid, rooted species-tree inference with 
    // SpeciesRax, and gene-tree species-tree reconciliation and estimation 
    // of gene family duplication transfer and loss with GeneRax. 
    
    // Do this for both the core and non-core gene families.
    // All outputs are needed for species tree inference, but not for the 
    // remainder. 
    SPECIES_TREE_PREP (
        ch_all_core_trees,
        ch_all_core_msas
    )
    .set { ch_core_spptree_prep }
    
    ch_core_treefile = ch_core_spptree_prep.treefile
    ch_core_families = ch_core_spptree_prep.families
    ch_core_generax_map = ch_core_spptree_prep.generax_map
    ch_asteroid_map = ch_core_spptree_prep.asteroid_map
    
    GENE_TREE_PREP (
        ch_all_rem_trees,
        ch_all_rem_msas
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
    ASTEROID (
        ch_core_treefile,
        ch_asteroid_map
    )
    .spp_tree
    .set { ch_asteroid }
    
    //
    // MODULE: SPECIESRAX
    // Now infer the rooted species tree with SpeciesRax,
    // reconcile gene family trees, and infer per-family 
    // rates of gene-family duplication, transfer, and loss
    //
    SPECIESRAX (
        ch_asteroid,
        ch_core_generax_map,
        ch_all_core_trees,
        ch_all_core_msas,
        ch_core_families
    )
    .speciesrax_tree
    .set { ch_speciesrax }
    
    // Run again, but this time only using the GeneRax component, 
    // reconciling gene family trees with the rooted species tree 
    // inferred from SpeciesRax for all remaining gene families
    GENERAX (
        ch_speciesrax,
        ch_rem_generax_map,
        ch_all_rem_trees,
        ch_all_rem_msas,
        ch_rem_families
    )
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//workflow.onComplete {
//    if (params.email || params.email_on_fail) {
//        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
//    }
//    NfcoreTemplate.summary(workflow, params, log)
//    if (params.hook_url) {
//       NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
//    }
//}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

