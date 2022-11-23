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
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Full samplesheet not specified!' }
if (params.test_data) { ch_test_data = file(params.test_data) } else { exit 1, 'Test samplesheet not specified!' }
if (params.fasta_dir) { ch_fa_dir = params.fasta_dir } else { exit 1, 'Fasta directory not specified!' }
if (params.test_fasta_dir) { ch_test_fa_dir = params.test_fasta_dir } else { exit 1, 'Test fasta directory not specified!' }
if (params.mcl_inflation) { ch_mcl_inflation = Channel.of(params.mcl_inflation) } else { exit 1, 'MCL Inflation parameter(s) not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
//ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


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
include { ORTHOFINDER_PREP                          } from '../modules/local/orthofinder_prep'
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_TEST } from '../modules/local/orthofinder_prep'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_TEST   } from '../modules/local/orthofinder_mcl'
include { ORTHOFINDER_MCL                           } from '../modules/local/orthofinder_mcl'
include { ANNOTATE_UNIPROT                          } from '../modules/local/annotate_uniprot'
include { COGEQC                                    } from '../modules/local/cogeqc'
include { SELECT_INFLATION                          } from '../modules/local/select_inflation'
include { CLIPKIT                                   } from '../modules/local/clipkit'
include { SPECIES_TREE_PREP                         } from '../modules/local/species_tree_prep'
include { ASTEROID                         } from '../modules/local/asteroid'
include { SPECIESRAX                         } from '../modules/local/speciesrax'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BUSCO                                 } from '../modules/nf-core/busco/main'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_TEST } from '../modules/nf-core/diamond/blastp/main'
include { DIAMOND_BLASTP                        } from '../modules/nf-core/diamond/blastp/main'
include { MAFFT                                 } from '../modules/nf-core/mafft/main'
include { IQTREE                                 } from '../modules/nf-core/iqtree/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//def meta = ch_spp
// This could be put into either the config file or specified via commandline. 

// Function to get list of [ meta, [ file ] ]
def create_og_channel(LinkedHashMap row) {
    // create meta map
    def meta  = [:]
        meta.og   = row.orthogroup
        meta.nspp = row.num_spp
        meta.copy_num = row.mean_copy_num
    // add path(s) of the OG file to the meta map
    def og_meta = []
        og_meta = [ meta, [ file(row.file) ] ] 
    return og_meta
}

workflow PHYLORTHOLOGY {
    
    ch_inflation = ch_mcl_inflation.toList().flatten()
    ch_versions = Channel.empty()
    //ch_fa = Channel.empty()
    //ch_dmd = Channel.empty()
    
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
    // First at the taxonomically shallow scale
    //if (params.busco_lineages_path) { 
    //    ch_busco_dat = file(params.busco_lineages_path)
    //    BUSCO (
    //        ch_busco_full,
    //        "taxon_specific",
    //        ch_busco_dat
    //    )
    //} else { 
    //    BUSCO (
    //        ch_busco_full,
    //        "taxon_specific",
    //        []
    //    )
    //}
    // Then at the broad scale (i.e. eukaryota)
    //if (params.busco_lineages_path) { 
    //    ch_busco_dat = file(params.busco_lineages_path)
    //    BUSCO (
    //        ch_busco_full,
    //        "eukaryota",
    //        ch_busco_dat
    //    )
    //} else { 
    //    BUSCO (
    //        ch_busco_full,
    //        "eukaryota",
    //       []
    //    )
    //}

    //
    // MODULE: All-v-All diamond/blastp
    //
    // Set the publishdir so that the similarity scores can be accessed by 
    // orthofinder in the next step (MCL clustering into orthogroups).
    ch_blastp_test = DIAMOND_BLASTP_TEST (
        ch_orthofinder_test,
        ch_test_dmd,
        "txt",
        "true",
        []
    )
    .txt
    
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

    //ch_similarities_remaining = ch_blastp.mix(ch_blastp).collect()
    //ch_similarities = ch_similarities_test.mix(ch_similarities_remaining)
    
    // First determine the optimal MCL inflation parameter, and then 
    // subsequently use this for full orthogroup inference. 
    ch_mcl = ORTHOFINDER_MCL_TEST (
        ch_inflation,
        ch_similarities_test,
        "true",
        [],
        []
    )
    .og_fpath
    
    //
    // MODULE: Run an R-script that applies cogqc to assess orthogroup inference
    //         accuracy/performance.
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
        "7",
        "1.5"
    )
    .core_ogs
    .splitCsv ( header:true, sep:',' )
    .map { create_og_channel(it) }
 
    ch_og_msas = MAFFT (
        ch_orthogroups
    )
    .msas
    
    ch_trimmed_msas = CLIPKIT (
        ch_og_msas
    )
    .trimmed_msas
    
    ch_gene_trees = IQTREE (
        ch_trimmed_msas,
        []
    )
    .phylogeny
    
    // Collect these gene family trees and alignments;
    // they will be used for unrooted species tree inference 
    // with Asteroid and downstream analysis with GeneRax and 
    // SpeciesRax
    // First trees....
    ch_gene_trees
    .branch {
        meta, phylogeny ->
            trees  : phylogeny
                return phylogeny
    }
    .collect()
    .set { ch_all_trees } 
    
    ch_trimmed_msas
    .branch {
        meta, trimmed_msas ->
            msas  : trimmed_msas
                return trimmed_msas
    }
    .collect()
    .set { ch_all_msas } 
    
    // Now, go ahead and prepare input files for initial unrooted species 
    // tree inference with Asteroid, rooted species-tree inference with 
    // SpeciesRax, and gene-tree species-tree reconciliation and estimation 
    // of gene family duplication transfer and loss with GeneRax. 
    SPECIES_TREE_PREP (
        ch_all_trees,
        ch_all_msas
    )
    .set { ch_spp_tree_prep }
    
    ch_treefile = ch_spp_tree_prep.treefile
    ch_families = ch_spp_tree_prep.families
    ch_generax_map = ch_spp_tree_prep.generax_map
    ch_asteroid_map = ch_spp_tree_prep.asteroid_map
    
    // Alrighty, now let's infer an intial, unrooted species tree using 
    // Asteroid
    ASTEROID (
        ch_treefile,
        ch_asteroid_map
    )
    .spp_tree
    .set { ch_asteroid }
    
    // Now infer the rooted species tree with SpeciesRax,
    // reconcile gene family trees, and infer per-family 
    // rates of gene-family duplication, transfer, and loss
    SPECIESRAX (
        ch_asteroid,
        ch_generax_map,
        ch_all_trees,
        ch_all_msas,
        ch_families
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

