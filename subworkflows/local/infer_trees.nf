include { PREQUAL } from '../../modules/local/prequal'

// Legacy single-aligner includes (used when adaptive_align = false)
if (!params.adaptive_align) {
    if (params.aligner == "witch") {
        include { WITCH as ALIGN_SEQS             } from '../../modules/local/witch'
    } else if (params.aligner == "famsa") {
        include { FAMSA as ALIGN_SEQS             } from '../../modules/local/famsa'
    } else {
        include { MAFFT as ALIGN_SEQS             } from '../../modules/nf-core-modified/mafft'
    }
}

// Adaptive alignment includes (used when adaptive_align = true)
if (params.adaptive_align) {
    include { MAFFT_ADAPTIVE                      } from '../../modules/local/mafft_adaptive'
    include { WITCH as WITCH_TIER2                } from '../../modules/local/witch'
    include { FAMSA as FAMSA_TIER3                } from '../../modules/local/famsa'
    include { FAMSA as FAMSA_FALLBACK             } from '../../modules/local/famsa'
}

if (params.msa_trimmer == "clipkit") {
    include { CLIPKIT as TRIM_MSAS            } from '../../modules/local/clipkit'
} else if (params.msa_trimmer == 'cialign') {
    include { CIALIGN as TRIM_MSAS            } from '../../modules/local/cialign'
}

if (params.tree_method == "iqtree") {
    include { IQTREE as TREES                 } from '../../modules/nf-core-modified/iqtree'
} else {
    include { FASTTREE as TREES               } from '../../modules/local/fasttree'
}

// FastTree fallback for when IQ-TREE fails (only in fallback mode)
if (params.iqtree_fasttree_fallback) {
    include { FASTTREE as FASTTREE_FALLBACK   } from '../../modules/local/fasttree'
}

workflow INFER_TREES {
    take:
    fas

    main:
    versions = Channel.empty()

    // Pre-alignment masking of non-homologous segments
    PREQUAL(fas)
    versions = versions.mix(PREQUAL.out.versions)

    if (params.adaptive_align) {
        // ── Size-adaptive alignment routing ──────────────────────────────
        // Count sequences per OG and branch by tier thresholds
        prequal_with_count = PREQUAL.out.masked
            .map { meta, fasta ->
                def n_seq = fasta.text.count('>')
                [meta + [n_seq: n_seq], fasta]
            }

        // Branch into tiers based on sequence count
        prequal_with_count.branch {
            tier1: it[0].n_seq <= params.align_tier1_max
            tier2: it[0].n_seq <= params.align_tier2_max
            tier3: true
        }.set { tiered }

        // Tier 1: MAFFT E-INS-i or L-INS-i (small families, ≤300 seqs)
        MAFFT_ADAPTIVE(tiered.tier1)
        versions = versions.mix(MAFFT_ADAPTIVE.out.versions)

        // Tier 2: WITCH (medium families, 301–3000 seqs)
        WITCH_TIER2(tiered.tier2)
        versions = versions.mix(WITCH_TIER2.out.versions)

        // Tier 3: FAMSA2 with accuracy flags (large families, >3000 seqs)
        FAMSA_TIER3(tiered.tier3)
        versions = versions.mix(FAMSA_TIER3.out.versions)

        // Detect tier 1 failures: inputs that didn't produce alignments
        tier1_failed = tiered.tier1
            .map { meta, fasta -> [meta.og, meta, fasta] }
            .join(
                MAFFT_ADAPTIVE.out.msas.map { meta, aln -> [meta.og, aln] },
                remainder: true
            )
            .filter { it[3] == null }  // No alignment = MAFFT failed
            .map { og, meta, fasta, aln -> [meta, fasta] }

        // Detect tier 2 failures: inputs that didn't produce alignments
        tier2_failed = tiered.tier2
            .map { meta, fasta -> [meta.og, meta, fasta] }
            .join(
                WITCH_TIER2.out.msas.map { meta, aln -> [meta.og, aln] },
                remainder: true
            )
            .filter { it[3] == null }  // No alignment = WITCH failed
            .map { og, meta, fasta, aln -> [meta, fasta] }

        // Run FAMSA fallback on all tier 1/2 failures
        FAMSA_FALLBACK(tier1_failed.mix(tier2_failed))
        versions = versions.mix(FAMSA_FALLBACK.out.versions)

        // Combine all alignment outputs
        all_msas = MAFFT_ADAPTIVE.out.msas
            .mix(WITCH_TIER2.out.msas)
            .mix(FAMSA_TIER3.out.msas)
            .mix(FAMSA_FALLBACK.out.msas)

        all_map_links = MAFFT_ADAPTIVE.out.map_link
            .mix(WITCH_TIER2.out.map_link)
            .mix(FAMSA_TIER3.out.map_link)
            .mix(FAMSA_FALLBACK.out.map_link)

    } else {
        // ── Legacy single-aligner path ──────────────────────────────────
        ALIGN_SEQS(PREQUAL.out.masked)
        versions = versions.mix(ALIGN_SEQS.out.versions)

        all_msas = ALIGN_SEQS.out.msas
        all_map_links = ALIGN_SEQS.out.map_link
    }

    if (params.msa_trimmer != 'none') {
        TRIM_MSAS(all_msas)
        versions = versions.mix(TRIM_MSAS.out.versions)
        map_link = TRIM_MSAS.out.map_link
        cleaned_msas = TRIM_MSAS.out.cleaned_msas
    } else {
        map_link = all_map_links
        cleaned_msas = all_msas
    }

    // Run primary tree inference (IQTREE in fallback mode, otherwise based on tree_method)
    TREES(cleaned_msas, params.tree_model)
    versions = versions.mix(TREES.out.versions)

    if (params.iqtree_fasttree_fallback && params.tree_method == "iqtree") {
        // Detect failed alignments by finding inputs that didn't produce trees
        failed_alignments = cleaned_msas
            .map { meta, aln -> [meta.og, meta, aln] }
            .join(
                TREES.out.phylogeny.map { meta, tree -> [meta.og, tree] },
                remainder: true
            )
            .filter { it[3] == null }  // No tree = IQ-TREE failed
            .map { og, meta, aln, tree -> [meta, aln] }

        // Run FastTree on failed alignments
        FASTTREE_FALLBACK(failed_alignments, params.tree_model)
        versions = versions.mix(FASTTREE_FALLBACK.out.versions)

        // Combine successful IQ-TREE trees with FastTree fallback trees
        phylogeny = TREES.out.phylogeny.mix(FASTTREE_FALLBACK.out.phylogeny)
    } else {
        phylogeny = TREES.out.phylogeny
    }

    emit:
    phylogeny
    map_link
    cleaned_msas
    versions
}
