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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
//ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
//ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
//ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK    } from '../subworkflows/local/input_check'

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
    CREATE CHANNELS: SPECIES NAME, FILE NAME, SHALLOW LINEAGE SPEC, SECOND LINEAGE SPEC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// This could probably go in the config file - also not certain it's being used. Seems to still be downloading data?
//buscoDatChannel = Channel.fromPath( '~/environment/resources/busco/busco_databases_v5.4.3/lineages/')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//def meta = ch_spp
// This could be put into either the config file or specified via commandline. 

workflow PHYLORTHOLOGY {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    ch_prots = INPUT_CHECK (
        ch_input
    )
    .prots

    ch_prots
    .map {
        meta, prots ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, prots ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, prots ->
            proteomes  : prots.size() == 1
                return [ meta, prots.flatten() ]
    }
    .set { ch_busco }

    ch_prots
    .map {
        meta, prots ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, prots ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, prots ->
            proteomes  : prots.size() == 1
                return prots.flatten()
    }
    .set { ch_diamond }
    
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


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
    if (params.busco_lineages_path) { 
        ch_busco_dat = file(params.busco_lineages_path)
        BUSCO (
            ch_busco,
            "eukaryota",
            ch_busco_dat
        )
    } else { 
        BUSCO (
            ch_busco,
            "eukaryota",
           []
        )
    }
    
    //
    // MODULE: Make diamond databases
    //
    // need to extract fasta file paths from ch_fasta for input into diamond. 
    ch_dbs_list = DIAMOND_MAKEDB (
        ch_diamond
    )
    .db
    
    ch_dbs = ch_dbs_list.mix(ch_dbs_list).collect()

    //
    // MODULE: All-v-All diamond/blastp
    //
    ch_blast = DIAMOND_BLASTP (
        ch_busco,
        ch_dbs,
        "txt",
        []
    )
    .txt
    
    
    

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

