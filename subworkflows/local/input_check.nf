//
// Check input samplesheet and get proteome channels
//

include { SAMPLESHEET_CHECK          } from '../../modules/local/samplesheet_check' params(params)

include { MCL_TEST_SAMPLESHEET_CHECK } from '../../modules/local/mcl_test_samplesheet_check' params(params)


workflow INPUT_CHECK {
    take:
    s3_dir
    complete_samplesheet // file: /path/to/complete_samplesheet.csv
    mcl_test_samplesheet // file: /path/to/mcl_test_samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( s3_dir, complete_samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_prots_channel(it) }
        .set { prots }
    
    MCL_TEST_SAMPLESHEET_CHECK ( s3_dir, mcl_test_samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_prots_channel(it) }
        .set { mcl_test_prots }
    
    emit:
    prots                                                  // channel: [ val(meta), [ prots ] ]
    mcl_test_prots                                         // channel: [ val(meta), [ mcl_test_prots ] ]
    all_data_prep = SAMPLESHEET_CHECK.out.of_prep          // path to a file that specifies where proteomes were downloaded into - to prep for orthofinder
    mcl_test_prep = MCL_TEST_SAMPLESHEET_CHECK.out.of_prep // path to a file that specifies where proteomes were downloaded into - to prep for orthofinder
    versions = SAMPLESHEET_CHECK.out.versions              // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ file ] ]
def create_prots_channel(LinkedHashMap row) {
    // create meta map
    def meta  = [:]
        meta.id   = row.species
        meta.taxon = row.taxonomy
        meta.shallow = row.shallow
        meta.broad = row.broad
        meta.mode = row.mode
        meta.uniprot = row.uniprot
        meta.mcl_test = row.mcl_test

    // add path(s) of the proteome file to the meta map
    def prots_meta = []
        prots_meta = [ meta, [ file(row.fpath) ] ] 
    return prots_meta
}