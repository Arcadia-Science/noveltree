if (params.aligner == "witch") {
    include { WITCH as ALIGN_SEQS             } from '../../modules/local/witch'
} else if (params.aligner == "famsa") {
    include { FAMSA as ALIGN_SEQS             } from '../../modules/local/famsa'
} else {
    include { MAFFT as ALIGN_SEQS             } from '../../modules/nf-core-modified/mafft'
}

if (params.msa_trimmer == "clipkit") {
    include { CLIPKIT as TRIM_MSAS            } from '../../modules/local/clipkit'
} else if (params.msa_trimmer == 'cialign') {
    include { CIALIGN as TRIM_MSAS            } from '../../modules/local/cialign'
}

// Tree inference: IQTREE when using iqtree method or fallback mode, otherwise FastTree
if (params.iqtree_fasttree_fallback || params.tree_method == "iqtree") {
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

    ALIGN_SEQS(fas)
    versions = versions.mix(ALIGN_SEQS.out.versions)

    if (params.msa_trimmer != 'none') {
        TRIM_MSAS(ALIGN_SEQS.out.msas)
        versions = versions.mix(TRIM_MSAS.out.versions)
        map_link = TRIM_MSAS.out.map_link
        cleaned_msas = TRIM_MSAS.out.cleaned_msas
    } else {
        map_link = ALIGN_SEQS.out.map_link
        cleaned_msas = ALIGN_SEQS.out.msas
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
