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
            remote: it[0].needs_download == true
            local:  true
        }
        .set { ch_branched }

    // Local file entries: [ val(meta), path(fasta) ]
    ch_local = ch_branched.local

    // Remote entries (URLs, NCBI accessions, UniProt IDs): [ val(meta), val(source_string) ]
    ch_remote = ch_branched.remote

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
    remote_prots   = ch_remote                   // channel: [ val(meta), val(source_string) ]
    mcl_test_prots                               // channel: [ val(meta), path(fasta) ]
    uniprot_prots                                // channel: [ val(meta), path(fasta) ]
    annotation_prots                             // channel: [ val(meta), path(fasta) ]
    complete_samplesheet = SAMPLESHEET_CHECK.out.csv
    versions = SAMPLESHEET_CHECK.out.versions    // channel: [ versions.yml ]
}

// Function to get list of [meta, file_or_source_string]
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

    // Detect source type from file column value
    def is_url = row.file.startsWith('http://') || row.file.startsWith('https://') ||
                 row.file.startsWith('ftp://') || row.file.startsWith('s3://')
    def is_ncbi = row.file ==~ /^GCF_\d+(\.\d+)?$/
    def is_uniprot = row.file ==~ /^UP\d{9,}$/

    if (is_ncbi) {
        meta.source_type = 'ncbi_refseq'
        // RefSeq proteomes always include isoforms — auto-override
        meta.isoform = 'yes'
    } else if (is_uniprot) {
        meta.source_type = 'uniprot'
    } else if (is_url) {
        meta.source_type = 'url'
    } else {
        meta.source_type = 'local'
    }
    meta.needs_download = (meta.source_type != 'local')

    if (meta.needs_download) {
        return [meta, row.file]
    } else {
        return [meta, file(row.file)]
    }
}
