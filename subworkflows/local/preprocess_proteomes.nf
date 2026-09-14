//
// Preprocess proteomes: optional TransDecoder, isoform filtering, and quality cleanup
//

include { TRANSDECODER        } from '../../modules/local/transdecoder'
include { FILTER_ISOFORMS     } from '../../modules/local/filter_isoforms'
include { FILTER_ISOFORMS_GFF3 } from '../../modules/local/filter_isoforms_gff3'
include { PREPROCESS_PROTEOME } from '../../modules/local/preprocess_proteome'

workflow PREPROCESS_PROTEOMES {
    take:
    ch_prots            // [ val(meta), path(fasta) ]
    min_protein_length  // val
    ch_annotations     // [ val(species), path(gff3) ]

    main:

    // Step 1: TransDecoder — only for species with meta.transdecoder == 'yes'
    ch_needs_transdecoder = ch_prots.filter { it[0].transdecoder == 'yes' }
    ch_skip_transdecoder  = ch_prots.filter { it[0].transdecoder != 'yes' }

    TRANSDECODER(ch_needs_transdecoder)

    ch_after_transdecoder = TRANSDECODER.out.translated.mix(ch_skip_transdecoder)

    // Step 2: Isoform filtering — for species with isoform == 'yes' OR transdecoder == 'yes'
    ch_needs_isofilter = ch_after_transdecoder.filter {
        it[0].filter_isoforms == 'yes' || it[0].transdecoder == 'yes'
    }
    ch_skip_isofilter = ch_after_transdecoder.filter {
        it[0].filter_isoforms != 'yes' && it[0].transdecoder != 'yes'
    }

    FILTER_ISOFORMS(ch_needs_isofilter.filter { meta, fasta -> !meta.has_gff3 })
    ch_gff_inputs = ch_needs_isofilter
        .filter { meta, fasta -> meta.has_gff3 }
        .map { meta, fasta -> tuple(meta.id, meta, fasta) }
        .join(ch_annotations, failOnDuplicate: true, failOnMismatch: true)
        .map { id, meta, fasta, gff3 -> tuple(meta, fasta, gff3) }
    FILTER_ISOFORMS_GFF3(ch_gff_inputs)

    ch_gene_filtered = FILTER_ISOFORMS_GFF3.out.filtered
        .map { meta, fasta -> tuple(meta + [gene_filtered: true], fasta) }
    ch_after_isofilter = FILTER_ISOFORMS.out.filtered.mix(ch_skip_isofilter, ch_gene_filtered)

    // Step 3: General preprocessing — always runs for every species
    // (stop codon cleanup, rare amino acids, length filter, CD-HIT)
    PREPROCESS_PROTEOME(ch_after_isofilter, min_protein_length)

    emit:
    preprocessed = PREPROCESS_PROTEOME.out.preprocessed  // [ val(meta), path(fasta) ]
    isoform_ledger = FILTER_ISOFORMS_GFF3.out.ledger
    isoform_qc = FILTER_ISOFORMS_GFF3.out.qc
}
