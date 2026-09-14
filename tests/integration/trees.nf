include { INFER_GENE_TREES } from '../../subworkflows/local/infer_gene_trees'

workflow {
    inputs = Channel.fromPath("${params.fixtures}/*.fa")
        .map { fasta ->
            tuple([og: fasta.baseName, n_seq: fasta.baseName == 'OGlarge' ? 100 : 4, max_len: 40],
                  fasta, file("${params.fixtures}/${fasta.baseName}.map.link"))
        }
    INFER_GENE_TREES(inputs)
    INFER_GENE_TREES.out.original_phylogeny.view { meta, tree -> "TREE ${meta.og}: ${tree.name}" }
}
