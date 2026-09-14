include { INPUT_CHECK } from '../../subworkflows/local/input_check'
include { PREPARE_INPUTS } from '../../subworkflows/local/prepare_inputs'
workflow {
    INPUT_CHECK(file(params.input))
    PREPARE_INPUTS(INPUT_CHECK.out.remote_prots, INPUT_CHECK.out.local_prots)
    PREPARE_INPUTS.out.renamed_prots.view { meta, fasta -> "PROTEOME ${meta.id}: ${fasta.name}" }
}
