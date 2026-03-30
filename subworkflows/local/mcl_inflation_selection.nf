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
            exit 1, 'Samplesheet must include samples with has_uniprot_ids=yes and include_in_mcl_test=yes when performing MCL parameter selection.'
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

        // Build species ID → name map for descriptive blast directory names
        ch_test_spp_id_map = ORTHOFINDER_PREP_TEST.out.sppIDs
            .map { file ->
                def nameMap = [:]
                file.text.trim().split('\n').each { line ->
                    def parts = line.split(': ')
                    nameMap[parts[0].trim()] = parts[1].replace('.fa', '')
                }
                nameMap
            }

        // Run for the test set (used to determine the best value of the MCL
        // inflation parameter)
        DIAMOND_BLASTP_TEST(
            mcl_test_prots,
            ORTHOFINDER_PREP_TEST.out.fastas.flatten(),
            ORTHOFINDER_PREP_TEST.out.diamonds.flatten(),
            "txt",
            "true",
            ch_test_spp_id_map
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
            "mcl_test_dataset",
            [],   // samplesheet — not used for mcl_test_dataset
            0,    // min_num_seqs — not used for mcl_test_dataset
            0,    // min_num_spp — not used for mcl_test_dataset
            0,    // min_prop_spp_for_spptree — not used for mcl_test_dataset
            0     // max_copy_num — not used for mcl_test_dataset
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

