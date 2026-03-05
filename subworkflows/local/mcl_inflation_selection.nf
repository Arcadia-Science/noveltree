include { ANNOTATE_UNIPROT                          } from '../../modules/local/annotate_uniprot'
include { COGEQC                                    } from '../../modules/local/cogeqc'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_TEST     } from '../../modules/nf-core-modified/diamond_blastp'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_TEST   } from '../../modules/local/orthofinder_mcl'
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_TEST } from '../../modules/local/orthofinder_prep'
include { SELECT_INFLATION                          } from '../../modules/local/select_inflation'

workflow MCL_INFLATION_SELECTION {
    take:
    mcl_test_prots
    annotation_prots
    mcl_inflation

    main:
        versions = Channel.empty()
        ch_inflation = Channel.fromList(mcl_inflation)
        mcl_test_prots_list = mcl_test_prots.collect { it[1] }

        annotation_prots.ifEmpty {
            exit 1, 'Samplesheet must include samples with UniProt annotations (uniprot column set to true) and mcl_test when performing MCL parameter selection.'
        }

        //
        // MODULE: Annotate UniProt Proteins
        //
        ANNOTATE_UNIPROT(annotation_prots)
            .cogeqc_annotations
            .collect()
            .set { ch_annotations }
        versions = versions.mix(ANNOTATE_UNIPROT.out.versions)

        ORTHOFINDER_PREP_TEST(mcl_test_prots_list, "mcl_test_dataset")

        // Run for the test set (used to determine the best value of the MCL
        // inflation parameter)
        DIAMOND_BLASTP_TEST(
            mcl_test_prots,
            ORTHOFINDER_PREP_TEST.out.fastas.flatten(),
            ORTHOFINDER_PREP_TEST.out.diamonds.flatten(),
            "txt",
            "true"
        )

        // First determine the optimal MCL inflation parameter, and then
        // subsequently use this for full orthogroup inference.
        ORTHOFINDER_MCL_TEST(
            ch_inflation,
            DIAMOND_BLASTP_TEST.out.txt.collect(),
            ORTHOFINDER_PREP_TEST.out.fastas,
            ORTHOFINDER_PREP_TEST.out.diamonds,
            ORTHOFINDER_PREP_TEST.out.sppIDs,
            ORTHOFINDER_PREP_TEST.out.seqIDs,
            "mcl_test_dataset"
        )

        COGEQC(
            ORTHOFINDER_MCL_TEST.out.inflation_dir,
            params.min_num_spp_per_og,
            ch_annotations
        )
        ch_cogeqc_summary = COGEQC.out.cogeqc_summary.collect()
        versions = versions.mix(COGEQC.out.versions)

        // Now, from these orthogroup summaries, select the best inflation parameter
        best_inflation = SELECT_INFLATION(ch_cogeqc_summary, params.min_num_spp_per_og)
            .best_inflation.text.trim().collect()
        versions = versions.mix(SELECT_INFLATION.out.versions)

    emit:
    best_inflation
    versions
}

