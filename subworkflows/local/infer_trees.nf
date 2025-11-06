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

if (params.tree_method == "iqtree") {
    include { IQTREE as TREES                 } from '../../modules/nf-core-modified/iqtree'
} else {
    include { FASTTREE as TREES               } from '../../modules/local/fasttree'
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

    TREES(cleaned_msas, params.tree_model)
    phylogeny = TREES.out.phylogeny
    versions = versions.mix(TREES.out.versions)

    emit:
    phylogeny
    map_link
    cleaned_msas
    versions
}
