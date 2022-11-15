//
// Check input samplesheet and get proteome channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_prots_channel(it) }
        .set { prots }

    emit:
    prots                                     // channel: [ val(meta), [ prots ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ file ] ]
def create_prots_channel(LinkedHashMap row) {
    // create meta map
    def meta  = [:]
    meta.id   = row.species
    meta.tax1 = row.tax1
    meta.tax2 = row.tax2
    meta.mode = row.mode
    meta.uniprot = row.uniprot

    // add path(s) of the proteome file to the meta map
    def prots_meta = []
    prots_meta = [ meta, [ file(row.file) ] ] // This is the part file(row.file) that's causing relative/absolute path problems. Seems to be interfacing with ch_input - specified by file(params.input) in a wierd way.
    // [ file(row.file) ] -> [/home/ubuntu/environment/github/~/environment/github/phylorthology/assets/nf-test/final-proteins/Lokiarchaeota_archaeon-test-proteome.fasta]]
    return prots_meta
}