//
// Orthogroup inference: MCL inflation selection, OrthoFinder prep, DIAMOND, MCL clustering
//

if (params.test_run_mcl) {
    include { MCL_INFLATION_SELECTION } from './mcl_inflation_selection'
}
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_ALL } from '../../modules/local/orthofinder_prep'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_ALL     } from '../../modules/nf-core-modified/diamond_blastp'
include { BUNDLE_BLAST_RESULTS as BUNDLE_BLAST_RESULTS_ALL } from '../../modules/local/bundle_blast_results'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_ALL   } from '../../modules/local/orthofinder_mcl'
include { CLEAN_ORTHOFINDER_CHECKPOINT } from '../../modules/local/orthofinder_checkpoint_cleanup'

workflow INFER_ORTHOGROUPS {
    take:
    renamed_prots        // [ val(meta), path(fasta) ]
    complete_prots_list  // value: collected list of fasta paths
    mcl_inflation        // list of inflation values
    complete_samplesheet // path: validated samplesheet

    main:
    //
    // MCL inflation parameter selection (opt-in via --test_run_mcl)
    //
    if (params.test_run_mcl) {
        // Ensure multiple inflation values are provided when testing
        if (mcl_inflation.size() < 2) {
            exit 1, '--test_run_mcl requires multiple --mcl_inflation values (e.g. --mcl_inflation "1.1,1.3,1.5,2.0,3.0")'
        }
        // Derive MCL test and annotation subsets from renamed files
        ch_renamed_mcl_test = renamed_prots.filter { it[0].include_in_mcl_test == 'yes' }
        ch_renamed_annotation = ch_renamed_mcl_test.filter { it[0].has_uniprot_ids == 'yes' }

        MCL_INFLATION_SELECTION(
            ch_renamed_mcl_test,
            ch_renamed_annotation,
            mcl_inflation
        )
        ch_best_inflation = MCL_INFLATION_SELECTION.out.best_inflation
    } else {
        ch_best_inflation = Channel.of(mcl_inflation.first())
    }

    //
    // Prepare directory structure and fasta files according to
    // OrthoFinder's preferred format for downstream MCL clustering
    //
    ORTHOFINDER_PREP_ALL(complete_prots_list, "complete_dataset")

    // Build species ID → name map from OrthoFinder's SpeciesIDs.txt
    ch_spp_id_map = ORTHOFINDER_PREP_ALL.out.sppIDs
        .map { file ->
            def nameMap = [:]
            file.text.trim().split('\n').each { line ->
                def parts = line.split(': ')
                nameMap[parts[0].trim()] = parts[1].replace('.fa', '')
            }
            nameMap
        }

    //
    // All-v-All diamond/blastp
    //
    DIAMOND_BLASTP_ALL(
        renamed_prots,
        ORTHOFINDER_PREP_ALL.out.fastas.flatten(),
        ORTHOFINDER_PREP_ALL.out.diamonds.flatten(),
        "txt",
        "false",
        ch_spp_id_map
    )

    // Bundle the all-v-all search results by query species before passing
    // them to OrthoFinder. Without this intermediate step, an N-species run
    // stages N^2 individual BLAST files into a single AWS Batch task. Large
    // datasets can overwhelm the generated Nextflow Bash staging wrapper.
    ch_blast_groups = DIAMOND_BLASTP_ALL.out.txt
        .map { blast ->
            def matcher = blast.name =~ /^Blast(\d+)_/
            if (!matcher.find()) {
                throw new IllegalArgumentException("Unexpected DIAMOND output name: ${blast.name}")
            }
            tuple("complete_${matcher.group(1)}", blast)
        }
        .groupTuple()

    BUNDLE_BLAST_RESULTS_ALL(ch_blast_groups)

    // Using the best-performing inflation parameter, infer orthogroups for
    // all samples. Also runs chimera detection and orthogroup filtering.
    ORTHOFINDER_MCL_ALL(
        ch_best_inflation,
        BUNDLE_BLAST_RESULTS_ALL.out.archive.collect(),
        ORTHOFINDER_PREP_ALL.out.fastas,
        ORTHOFINDER_PREP_ALL.out.diamonds,
        ORTHOFINDER_PREP_ALL.out.sppIDs,
        ORTHOFINDER_PREP_ALL.out.seqIDs,
        "complete_dataset",
        complete_samplesheet,
        params.min_num_seq_per_og,
        params.min_num_spp_per_og,
        params.min_prop_spp_for_spptree,
        params.max_copy_num_spp_tree
    )

    // This receipt is emitted only after Nextflow has successfully finalized
    // and stored every ORTHOFINDER_MCL output. The cleanup task therefore
    // stages one tiny file and cannot erase recovery data during finalization.
    CLEAN_ORTHOFINDER_CHECKPOINT(ORTHOFINDER_MCL_ALL.out.checkpoint_receipt)

    // Parse one compact manifest and join its resource metadata to the FASTA
    // paths by orthogroup. This avoids controller-side reads of every FASTA.
    ch_og_metadata = ORTHOFINDER_MCL_ALL.out.og_metadata
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                og: row.orthogroup,
                n_seq: row.n_seq as int,
                max_len: row.max_len as int,
                n_species: row.n_species as int,
            ]
            tuple(row.orthogroup, row.family_set, meta)
        }

    ch_spptree_metadata = ch_og_metadata
        .filter { og, family_set, meta -> family_set == 'species_tree' }
        .map { og, family_set, meta -> tuple(og, meta) }
    ch_genetree_metadata = ch_og_metadata
        .filter { og, family_set, meta -> family_set == 'gene_tree' }
        .map { og, family_set, meta -> tuple(og, meta) }

    ch_spptree_fas = ORTHOFINDER_MCL_ALL.out.spptree_fas
        .flatten()
        .map { fasta -> tuple(fasta.simpleName, fasta) }
        .join(ch_spptree_metadata)
        .map { og, fasta, meta -> tuple(meta, fasta) }
    ch_genetree_fas = ORTHOFINDER_MCL_ALL.out.genetree_fas
        .flatten()
        .map { fasta -> tuple(fasta.simpleName, fasta) }
        .join(ch_genetree_metadata)
        .map { og, fasta, meta -> tuple(meta, fasta) }

    // --test mode: keep only gene families that contain ALL species in the dataset.
    // This dramatically reduces the number of OGs for quick end-to-end smoke tests.
    if (params.test_run) {
        ch_total_species = renamed_prots.collect { it[0].id }
            .map { ids -> ids.unique().size() }

        ch_spptree_fas = ch_spptree_fas
            .combine(ch_total_species)
            .filter { meta, fasta, n_spp ->
                meta.n_species >= n_spp
            }
            .map { meta, fasta, n_spp -> [meta, fasta] }

        ch_genetree_fas = ch_genetree_fas
            .combine(ch_total_species)
            .filter { meta, fasta, n_spp ->
                meta.n_species >= n_spp
            }
            .map { meta, fasta, n_spp -> [meta, fasta] }
    }

    emit:
    spptree_fas   = ch_spptree_fas
    genetree_fas  = ch_genetree_fas
    inflation_dir = ORTHOFINDER_MCL_ALL.out.inflation_dir
}
