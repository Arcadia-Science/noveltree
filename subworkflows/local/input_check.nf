//
// Check input samplesheet and get proteome channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    complete_samplesheet // file: /path/to/complete_samplesheet.csv
    data_dir
    
    main:
    SAMPLESHEET_CHECK ( complete_samplesheet, data_dir )
    
    SAMPLESHEET_CHECK.out
        .complete_csv
        .splitCsv ( header:true, sep:',' )
        .map { create_prots_channel(it) }
        .set { complete_prots }
        
    SAMPLESHEET_CHECK.out
        .mcl_test_csv
        .splitCsv ( header:true, sep:',' )
        .map { create_prots_channel(it) }
        .set { mcl_test_prots }
    
    emit:
    complete_prots                            // channel: [ val(meta), [ complete_prots ] ]
    mcl_test_prots                            // channel: [ val(meta), [ mcl_test_prots ] ]
    complete_samplesheet = SAMPLESHEET_CHECK.out.complete_csv
    complete_fastadir = SAMPLESHEET_CHECK.out.complete_fastadir
    mcl_test_fastadir = SAMPLESHEET_CHECK.out.mcl_test_fastadir
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
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
