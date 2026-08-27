//
// Download, preprocess, and rename proteome FASTA files
//

include { DOWNLOAD_INPUT } from '../../modules/local/download_input'
include { RENAME_FASTAS  } from '../../modules/local/rename_fastas'

if (params.preprocess) {
    include { PREPROCESS_PROTEOMES } from './preprocess_proteomes'
}

workflow PREPARE_INPUTS {
    take:
    remote_prots    // [ val(meta), val(source_string) ]
    local_prots     // [ val(meta), path(fasta) ]

    main:
    // Download any remote inputs (URLs, NCBI accessions, UniProt IDs), then merge with local files
    DOWNLOAD_INPUT(remote_prots)
    ch_complete_prots = local_prots.mix(DOWNLOAD_INPUT.out.downloaded)

    // Optional proteome preprocessing (TransDecoder, isoform filtering, quality cleanup)
    if (params.preprocess) {
        PREPROCESS_PROTEOMES(ch_complete_prots, params.min_protein_length)
        ch_to_rename = PREPROCESS_PROTEOMES.out.preprocessed
    } else {
        ch_to_rename = ch_complete_prots
    }

    // Rename FASTA files so filenames match normalized species names.
    // OrthoFinder uses filenames as species identifiers, so this ensures
    // all downstream tip labels use hyphens (e.g. Homo-sapiens_ProteinID).
    RENAME_FASTAS(ch_to_rename)

    emit:
    renamed_prots = RENAME_FASTAS.out.renamed
    protein_maps = RENAME_FASTAS.out.protein_map
    normalization_qc = RENAME_FASTAS.out.normalization_qc
}
