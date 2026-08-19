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
    include { FASTTREE as FASTTREE_TIER2      } from '../../modules/local/fasttree'
} else {
    include { FASTTREE as TREES               } from '../../modules/local/fasttree'
}

// FastTree fallback for when IQ-TREE fails (only in fallback mode)
if (params.iqtree_fasttree_fallback && params.tree_method == "iqtree") {
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

        // ── storeDir-aware pre-filter for alignment fallback ─────────────
        // On resume, OGs that previously fell back to FAMSA have only a
        // _famsa.fa in storeDir — no MAFFT/WITCH output.  Without this
        // filter, MAFFT/WITCH would rerun them (and fail again) before
        // routing to FAMSA_FALLBACK.  We check storeDir upfront.
        def alnStore = "${params.outdir}/alignments/original"

        tiered.tier1.branch { meta, fasta ->
            has_famsa: file("${alnStore}/${fasta.baseName}_famsa.fa").exists()
            needs_aligner: true
        }.set { tier1_routed }

        tiered.tier2.branch { meta, fasta ->
            has_famsa: file("${alnStore}/${fasta.baseName}_famsa.fa").exists()
            needs_aligner: true
        }.set { tier2_routed }

        // Tier 1: MAFFT E-INS-i or L-INS-i (small families, ≤200 seqs)
        MAFFT_TIER1(tier1_routed.needs_aligner)

        // Tier 2: WITCH (medium families, 201–1000 seqs by default)
        WITCH_TIER2(tier2_routed.needs_aligner)

        // Tier 3: FAMSA2 with accuracy flags (large families, >1000 seqs)
        FAMSA_TIER3(tiered.tier3)

        // Detect vanished OGs: present in input but absent from output.
        // Any failure (tool or infrastructure) is ignored by Nextflow,
        // producing no output for that OG — route it to FAMSA_FALLBACK.
        tier1_produced = MAFFT_TIER1.out.msas.map { meta, aln -> [meta.og, true] }
        newly_failed_tier1 = tier1_routed.needs_aligner
            .map { meta, fasta -> [meta.og, true] }
            .join(tier1_produced, remainder: true)
            .filter { it[2] == null }          // no match in output → vanished
            .map { it[0] }                     // OG id
            .join(tier1_routed.needs_aligner.map { meta, fasta -> [meta.og, meta, fasta] })
            .map { og, meta, fasta -> [meta, fasta] }

        tier2_produced = WITCH_TIER2.out.msas.map { meta, aln -> [meta.og, true] }
        newly_failed_tier2 = tier2_routed.needs_aligner
            .map { meta, fasta -> [meta.og, true] }
            .join(tier2_produced, remainder: true)
            .filter { it[2] == null }
            .map { it[0] }
            .join(tier2_routed.needs_aligner.map { meta, fasta -> [meta.og, meta, fasta] })
            .map { og, meta, fasta -> [meta, fasta] }

        // Run FAMSA fallback: fresh failures + pre-existing fallback OGs
        FAMSA_FALLBACK(
            newly_failed_tier1
                .mix(newly_failed_tier2)
                .mix(tier1_routed.has_famsa)
                .mix(tier2_routed.has_famsa)
        )

        // Combine successful primary alignments with fallback
        all_msas = MAFFT_TIER1.out.msas
            .mix(WITCH_TIER2.out.msas)
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

    if (params.tree_method == "iqtree") {
        // Direct very large families to FastTree instead of first waiting for
        // an expensive IQ-TREE failure or timeout.
        cleaned_msas.branch { meta, aln ->
            iqtree: (meta.n_seq as int) <= (params.tree_iqtree_max as int)
            fasttree: true
        }.set { tree_tiered }

        FASTTREE_TIER2(tree_tiered.fasttree, params.tree_model)

        if (params.iqtree_fasttree_fallback) {
            // ── storeDir-aware pre-filter ──────────────────────────────
            // On resume, OGs that previously fell back to FastTree have only
            // a _ft.newick in storeDir — no _iqt.newick. Send those OGs
            // straight to FASTTREE_FALLBACK, which hits storeDir immediately.
            def store = "${params.outdir}/gene_family_trees/original"
            tree_tiered.iqtree.branch { meta, aln ->
                has_ft: file("${store}/${aln.baseName}_ft.newick").exists()
                needs_iqtree: true
            }.set { routed }

            // Run IQ-TREE only on OGs that do not already have a FastTree result.
            TREES(routed.needs_iqtree, params.tree_model)

            // Detect vanished OGs: present in IQ-TREE input but absent from output.
            tree_produced = TREES.out.phylogeny.map { meta, tree -> [meta.og, true] }
            newly_failed = routed.needs_iqtree
                .map { meta, aln -> [meta.og, true] }
                .join(tree_produced, remainder: true)
                .filter { it[2] == null }
                .map { it[0] }
                .join(routed.needs_iqtree.map { meta, aln -> [meta.og, meta, aln] })
                .map { og, meta, aln -> [meta, aln] }

            FASTTREE_FALLBACK(newly_failed.mix(routed.has_ft), params.tree_model)

            phylogeny = TREES.out.phylogeny
                .mix(FASTTREE_FALLBACK.out.phylogeny)
                .mix(FASTTREE_TIER2.out.phylogeny)
        } else {
            TREES(tree_tiered.iqtree, params.tree_model)
            phylogeny = TREES.out.phylogeny.mix(FASTTREE_TIER2.out.phylogeny)
        }
    } else {
        // FastTree selected globally: run it on every family.
        TREES(cleaned_msas, params.tree_model)
        phylogeny = TREES.out.phylogeny
    }

    emit:
    phylogeny
    map_link
    cleaned_msas
}
