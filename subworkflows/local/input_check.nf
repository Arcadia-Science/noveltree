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

    emit:
    local_prots    = ch_local                    // channel: [ val(meta), path(fasta) ]
    remote_prots   = ch_remote                   // channel: [ val(meta), val(source_string) ]
    complete_samplesheet = SAMPLESHEET_CHECK.out.csv
    versions = SAMPLESHEET_CHECK.out.versions    // channel: [ versions.yml ]
}

// Function to get list of [meta, file_or_source_string]
def create_prots_channel(LinkedHashMap row) {
    // create meta map
    def meta  = [:]
        meta.id   = row.species
        meta.input_type = row.input_type
        meta.has_uniprot_ids = row.has_uniprot_ids
        meta.include_in_mcl_test = row.include_in_mcl_test
        meta.transdecoder = row.transdecoder
        meta.filter_isoforms = row.filter_isoforms
        meta.reference_proteome = row.reference_proteome
        meta.busco_shallow = row.busco_shallow
        meta.busco_broad = row.busco_broad

    // Detect source type from input_data column value
    def is_url = row.input_data.startsWith('http://') || row.input_data.startsWith('https://') ||
                 row.input_data.startsWith('ftp://') || row.input_data.startsWith('s3://')
    def is_ncbi = row.input_data ==~ /^GC[AF]_\d+(\.\d+)?$/
    def is_uniprot = row.input_data ==~ /^UP\d{9,}$/

    if (is_ncbi) {
        meta.source_type = 'ncbi_refseq'
        // RefSeq proteomes always include isoforms — auto-override
        meta.filter_isoforms = 'yes'
    } else if (is_uniprot) {
        meta.source_type = 'uniprot'
    } else if (is_url) {
        meta.source_type = 'url'
    } else {
        meta.source_type = 'local'
    }
    meta.needs_download = (meta.source_type != 'local')

    // Guard: TransDecoder is for nucleotide inputs only. If the input mode
    // is 'proteins', running ORF prediction will produce garbage. Auto-fix
    // and warn so the pipeline doesn't silently lose an entire species.
    if (meta.transdecoder == 'yes' && meta.input_type == 'proteins') {
        log.warn "Samplesheet has transdecoder=yes for '${meta.id}' but mode=proteins. " +
                 "TransDecoder requires nucleotide input — overriding to transdecoder=no."
        meta.transdecoder = 'no'
    }

    if (meta.needs_download) {
        return [meta, row.input_data]
    } else {
        return [meta, file(row.input_data)]
    }
}
