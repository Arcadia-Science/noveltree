//
// Check input samplesheet and get proteome channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    complete_samplesheet // file: /path/to/complete_samplesheet.csv

    main:
    SAMPLESHEET_CHECK(complete_samplesheet)

    SAMPLESHEET_CHECK.out
        .csv
        .splitCsv (header:true, sep:',')
        .map { create_prots_channel(it) }
        .set { complete_prots }

    complete_prots.filter {
        it[0].mcl_test == 'true'
    }.set { mcl_test_prots }

    emit:
    complete_prots                            // channel: [ val(meta), [ complete_prots ] ]
    mcl_test_prots                            // channel: [ val(meta), [ mcl_test_prots ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [meta, [file]]
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
        prots_meta = [meta, [file(row.file)]]
    return prots_meta
}
