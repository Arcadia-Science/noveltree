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
        .branch {
            urls:  it[0].is_url == true
            local: true
        }
        .set { ch_branched }

    // Local file entries: [ val(meta), path(fasta) ]
    ch_local = ch_branched.local

    // URL entries: [ val(meta), val(url_string) ]
    ch_urls = ch_branched.urls

    // Subset channels (from local files only — URL species join after download)
    ch_local.filter {
        it[0].mcl_test == 'true'
    }.set { mcl_test_prots }

    ch_local.filter {
        it[0].uniprot == 'true'
    }.set { uniprot_prots }

    mcl_test_prots.filter {
        it[0].uniprot == 'true'
    }.set { annotation_prots }

    emit:
    local_prots    = ch_local                    // channel: [ val(meta), path(fasta) ]
    url_prots      = ch_urls                     // channel: [ val(meta), val(url_string) ]
    mcl_test_prots                               // channel: [ val(meta), path(fasta) ]
    uniprot_prots                                // channel: [ val(meta), path(fasta) ]
    annotation_prots                             // channel: [ val(meta), path(fasta) ]
    complete_samplesheet = SAMPLESHEET_CHECK.out.csv
    versions = SAMPLESHEET_CHECK.out.versions    // channel: [ versions.yml ]
}

// Function to get list of [meta, file_or_url]
def create_prots_channel(LinkedHashMap row) {
    // create meta map
    def meta  = [:]
        meta.id   = row.species
        meta.taxon = row.taxonomy
        meta.shallow_db = row.shallow_db
        meta.broad_db = row.broad_db
        meta.mode = row.mode
        meta.uniprot = row.uniprot
        meta.mcl_test = row.mcl_test
        meta.annotate = (row.uniprot == "true") && (row.mcl_test == "true")
        meta.transdecoder = row.transdecoder ?: 'no'
        meta.isoform = row.isoform ?: 'no'
        meta.reference = row.reference ?: 'no'

    // Detect URLs — don't call file() on them to avoid Nextflow auto-staging
    def is_url = row.file.startsWith('http://') || row.file.startsWith('https://') ||
                 row.file.startsWith('ftp://') || row.file.startsWith('s3://')
    meta.is_url = is_url

    if (is_url) {
        return [meta, row.file]
    } else {
        return [meta, file(row.file)]
    }
}
