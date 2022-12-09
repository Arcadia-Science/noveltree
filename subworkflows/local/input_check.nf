//
// Check input samplesheet and get proteome channels
//

include { SAMPLESHEET_CHECK          } from '../../modules/local/samplesheet_check' //params(params)

workflow INPUT_CHECK {
    take:
    complete_samplesheet // file: /path/to/complete_samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( complete_samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_prots_channel(it) }
        .set { prots }
    
    emit:
    prots                                                  // channel: [ val(meta), [ prots ] ]
    all_data_prep = SAMPLESHEET_CHECK.out.of_prep          // path to a file that specifies where proteomes were downloaded into - to prep for orthofinder
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
        prots_meta = [ meta, [ file(row.file) ] ] 
    return prots_meta
}