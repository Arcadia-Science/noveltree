// Legacy single-aligner includes (used when aligner != 'adaptive')
if (params.aligner != 'adaptive') {
    if (params.aligner == "witch") {
        include { WITCH as ALIGN_SEQS             } from '../../modules/local/witch'
    } else if (params.aligner == "famsa") {
        include { FAMSA as ALIGN_SEQS             } from '../../modules/local/famsa'
    } else {
        include { MAFFT as ALIGN_SEQS             } from '../../modules/nf-core-modified/mafft'
    }
}

// Adaptive alignment includes (used when aligner == 'adaptive')
if (params.aligner == 'adaptive') {
    include { MAFFT_TIER1                         } from '../../modules/local/mafft_tier1'
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

workflow INFER_GENE_TREES {
    take:
    fas

    main:

    if (params.aligner == 'adaptive') {
        // ── Size-adaptive alignment routing ──────────────────────────────
        // Branch into tiers based on sequence count (meta.n_seq set in main.nf)
        fas.branch {
            tier1: it[0].n_seq <= params.align_tier1_max
            tier2: it[0].n_seq <= params.align_tier2_max
            tier3: true
        }.set { tiered }

        // Tier 1: MAFFT E-INS-i or L-INS-i (small families, ≤300 seqs)
        MAFFT_TIER1(tiered.tier1)

        // Tier 2: WITCH (medium families, 301–3000 seqs)
        WITCH_TIER2(tiered.tier2)

        // Tier 3: FAMSA2 with accuracy flags (large families, >3000 seqs)
        FAMSA_TIER3(tiered.tier3)

        // Filter real alignments (non-empty) from sentinels (0-byte)
        tier1_ok = MAFFT_TIER1.out.msas.filter { meta, aln -> aln.size() > 0 }
        tier2_ok = WITCH_TIER2.out.msas.filter { meta, aln -> aln.size() > 0 }

        // Detect failures: 0-byte sentinels, rejoin with input FASTA for fallback
        tier1_failed = MAFFT_TIER1.out.msas
            .filter { meta, aln -> aln.size() == 0 }
            .map { meta, aln -> [meta.og, meta] }
            .join(tiered.tier1.map { meta, fasta -> [meta.og, fasta] })
            .map { og, meta, fasta -> [meta, fasta] }

        tier2_failed = WITCH_TIER2.out.msas
            .filter { meta, aln -> aln.size() == 0 }
            .map { meta, aln -> [meta.og, meta] }
            .join(tiered.tier2.map { meta, fasta -> [meta.og, fasta] })
            .map { og, meta, fasta -> [meta, fasta] }

        // Run FAMSA fallback on all tier 1/2 failures
        FAMSA_FALLBACK(tier1_failed.mix(tier2_failed))

        // Combine successful primary alignments with fallback
        all_msas = tier1_ok
            .mix(tier2_ok)
            .mix(FAMSA_TIER3.out.msas)
            .mix(FAMSA_FALLBACK.out.msas)

        all_map_links = MAFFT_TIER1.out.map_link
            .mix(WITCH_TIER2.out.map_link)
            .mix(FAMSA_TIER3.out.map_link)
            .mix(FAMSA_FALLBACK.out.map_link)

    } else {
        // ── Legacy single-aligner path ──────────────────────────────────
        ALIGN_SEQS(fas)

        all_msas = ALIGN_SEQS.out.msas
        all_map_links = ALIGN_SEQS.out.map_link
    }

    if (params.msa_trimmer != 'none') {
        TRIM_MSAS(all_msas)
        // Filter out empty files (QC-failed OGs produce empty placeholders for storeDir)
        map_link = TRIM_MSAS.out.map_link.filter { meta, f -> f.size() > 0 }
        cleaned_msas = TRIM_MSAS.out.cleaned_msas.filter { meta, f -> f.size() > 0 }
    } else {
        map_link = all_map_links
        cleaned_msas = all_msas
    }

    // Run primary tree inference (IQTREE in fallback mode, otherwise based on tree_method)
    TREES(cleaned_msas, params.tree_model)

    if (params.iqtree_fasttree_fallback && params.tree_method == "iqtree") {
        // Filter real trees (non-empty) from sentinels (0-byte)
        real_trees = TREES.out.phylogeny.filter { meta, tree -> tree.size() > 0 }

        // Detect failures: 0-byte sentinels, rejoin with input alignment for fallback
        failed_alignments = TREES.out.phylogeny
            .filter { meta, tree -> tree.size() == 0 }
            .map { meta, tree -> [meta.og, meta] }
            .join(cleaned_msas.map { meta, aln -> [meta.og, aln] })
            .map { og, meta, aln -> [meta, aln] }

        // Run FastTree on failed alignments
        FASTTREE_FALLBACK(failed_alignments, params.tree_model)

        // Combine successful IQ-TREE trees with FastTree fallback trees
        phylogeny = real_trees.mix(FASTTREE_FALLBACK.out.phylogeny)
    } else {
        phylogeny = TREES.out.phylogeny
    }

    emit:
    phylogeny
    map_link
    cleaned_msas
}
