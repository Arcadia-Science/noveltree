include { WITCH } from '../../modules/local/witch'
workflow {
    WITCH(Channel.of(tuple([og: 'OGwitch', n_seq: 8, max_len: 80],
                           file("${params.fixtures}/OGwitch.fa"),
                           file("${params.fixtures}/OGwitch.map.link"))))
    WITCH.out.outcome.view { meta, outcome -> "WITCH OUTCOME: ${outcome.name}" }
}
