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
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta_dir) { ch_fa_dir = params.fasta_dir } else { exit 1, 'Fasta directory not specified!' }
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
include { INPUT_CHECK    } from '../subworkflows/local/input_check'

//
// MODULE
//
include { ORTHOFINDER_PREP           } from '../modules/local/orthofinder_prep'
include { ORTHOFINDER_MCL            } from '../modules/local/orthofinder_mcl'
include { ANNOTATE_UNIPROT           } from '../modules/local/annotate_uniprot'
include { COGEQC                     } from '../modules/local/cogeqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BUSCO                      } from '../modules/nf-core/busco/main'
include { DIAMOND_MAKEDB             } from '../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP             } from '../modules/nf-core/diamond/blastp/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//def meta = ch_spp
// This could be put into either the config file or specified via commandline. 

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
        ch_fa_dir
    )
    
    // Fasta files should be redirected into a channel set of filepaths emitted
    // separately, whereas the diamond databases for each species can be put
    // into a directory as they are now (as a comma seperated list emitted 
    // together).
    ch_fa = ORTHOFINDER_PREP.out.fa.flatten()
    ch_dmd = ORTHOFINDER_PREP.out.dmd.flatten()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_prots = INPUT_CHECK (
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


    ch_prots
    .branch {
        meta, prots ->
            proteomes  : prots.size() == 1
                return [ meta, prots.flatten() ]
    }
    .set { ch_busco }
    
    ch_of_fa = ch_busco.merge(ch_fa)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    //
    // MODULE: Annotate UniProt Proteins
    //
    ch_annotations = ANNOTATE_UNIPROT (
            ch_busco
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
    //        ch_busco,
    //        "taxon_specific",
    //        ch_busco_dat
    //    )
    //} else { 
    //    BUSCO (
    //        ch_busco,
    //        "taxon_specific",
    //        []
    //    )
    //}
    // Then at the broad scale (i.e. eukaryota)
    //if (params.busco_lineages_path) { 
    //    ch_busco_dat = file(params.busco_lineages_path)
    //    BUSCO (
    //        ch_busco,
    //        "eukaryota",
    //        ch_busco_dat
    //    )
    //} else { 
    //    BUSCO (
    //        ch_busco,
    //        "eukaryota",
    //       []
    //    )
    //}
    
    //
    // MODULE: Make diamond databases
    //
    // need to extract fasta file paths from ch_fasta for input into diamond. 
    //ch_dbs_list = DIAMOND_MAKEDB (
    //    ch_diamond
    //)
    //.db
    
    //ch_dbs = ch_dbs_list.mix(ch_dbs_list).collect()

    //
    // MODULE: All-v-All diamond/blastp
    //
    // Set the publishdir so that the similarity scores can be accessed by 
    // orthofinder in the next step (MCL clustering into orthogroups).
    ch_blastp = DIAMOND_BLASTP (
        ch_of_fa,
        ch_dmd,
        "txt",
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
    ch_similarities = ch_blastp.mix(ch_blastp).collect()

    ch_mcl = ORTHOFINDER_MCL (
        ch_inflation,
        ch_similarities
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

