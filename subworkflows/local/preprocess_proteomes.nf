//
// Preprocess proteomes: optional TransDecoder, isoform filtering, and quality cleanup
//

include { TRANSDECODER        } from '../../modules/local/transdecoder'
include { FILTER_ISOFORMS     } from '../../modules/local/filter_isoforms'
include { PREPROCESS_PROTEOME } from '../../modules/local/preprocess_proteome'

workflow PREPROCESS_PROTEOMES {
    take:
    ch_prots            // [ val(meta), path(fasta) ]
    min_protein_length  // val

    main:
    versions = Channel.empty()

    // Step 1: TransDecoder — only for species with meta.transdecoder == 'yes'
    ch_needs_transdecoder = ch_prots.filter { it[0].transdecoder == 'yes' }
    ch_skip_transdecoder  = ch_prots.filter { it[0].transdecoder != 'yes' }

    TRANSDECODER(ch_needs_transdecoder)
    versions = versions.mix(TRANSDECODER.out.versions)

    ch_after_transdecoder = TRANSDECODER.out.translated.mix(ch_skip_transdecoder)

    // Step 2: Isoform filtering — for species with isoform == 'yes' OR transdecoder == 'yes'
    ch_needs_isofilter = ch_after_transdecoder.filter {
        it[0].isoform == 'yes' || it[0].transdecoder == 'yes'
    }
    ch_skip_isofilter = ch_after_transdecoder.filter {
        it[0].isoform != 'yes' && it[0].transdecoder != 'yes'
    }

    FILTER_ISOFORMS(ch_needs_isofilter)
    versions = versions.mix(FILTER_ISOFORMS.out.versions)

    ch_after_isofilter = FILTER_ISOFORMS.out.filtered.mix(ch_skip_isofilter)

    // Step 3: General preprocessing — always runs for every species
    // (stop codon cleanup, rare amino acids, length filter, CD-HIT)
    PREPROCESS_PROTEOME(ch_after_isofilter, min_protein_length)
    versions = versions.mix(PREPROCESS_PROTEOME.out.versions)

    emit:
    preprocessed = PREPROCESS_PROTEOME.out.preprocessed  // [ val(meta), path(fasta) ]
    versions
}
