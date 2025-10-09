include { CIALIGN as TRIM_MSAS } from '../../modules/local/cialign'
include { FASTTREE             } from '../../modules/local/fasttree'
include { WITCH as ALIGN_SEQS  } from '../../modules/local/witch'

workflow INFER_TREES {
    take:
    fas

    main:
    versions = Channel.empty()

    ALIGN_SEQS(fas)
    versions = versions.mix(ALIGN_SEQS.out.versions)

    TRIM_MSAS(ALIGN_SEQS.out.cleaned_msas)
    versions = versions.mix(TRIM_MSAS.out.versions)
    map_link = TRIM_MSAS.out.map_link
    cleaned_msas = TRIM_MSAS.out.cleaned_msas

    // Pass dummy model parameter to match noveltree FASTTREE module signature
    // (FASTTREE doesn't use it, but needs it for interface compatibility with IQTREE)
    FASTTREE(cleaned_msas, 'LG+G4')
    phylogeny = FASTTREE.out.phylogeny
    versions = versions.mix(FASTTREE.out.versions)

    emit:
    phylogeny
    map_link
    cleaned_msas
    versions
}
