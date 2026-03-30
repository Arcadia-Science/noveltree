//
// Orthogroup inference: MCL inflation selection, OrthoFinder prep, DIAMOND, MCL clustering
//

if (params.test_run_mcl) {
    include { MCL_INFLATION_SELECTION } from './mcl_inflation_selection'
}
include { ORTHOFINDER_PREP as ORTHOFINDER_PREP_ALL } from '../../modules/local/orthofinder_prep'
include { DIAMOND_BLASTP as DIAMOND_BLASTP_ALL     } from '../../modules/nf-core-modified/diamond_blastp'
include { ORTHOFINDER_MCL as ORTHOFINDER_MCL_ALL   } from '../../modules/local/orthofinder_mcl'

// Function to get list of [meta, [file]]
def create_og_channel(Object inputs) {
    if (inputs instanceof String) {
        inputs = [inputs]
    }
    def metaList = []
    inputs.each { input ->
        def meta = [:]
        meta.og = input
        metaList.add(meta)
    }
    return metaList
}

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

    // Using the best-performing inflation parameter, infer orthogroups for
    // all samples. Also runs chimera detection and orthogroup filtering.
    ORTHOFINDER_MCL_ALL(
        ch_best_inflation,
        DIAMOND_BLASTP_ALL.out.txt.collect(),
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

    // Create meta maps for the two sets by providing the simple name of each orthogroup
    spptree_og_names = ORTHOFINDER_MCL_ALL.out.spptree_fas.map { file -> file.simpleName }
    spptree_og_map = spptree_og_names.map { create_og_channel(it) }.flatten()
    genetree_og_names = ORTHOFINDER_MCL_ALL.out.genetree_fas.map { file -> file.simpleName }
    genetree_og_map = genetree_og_names.map { create_og_channel(it) }.flatten()

    // Create the tuple of output fastas paired with the meta map,
    // injecting n_seq and max_len for downstream dynamic resource allocation.
    // max_len = longest sequence in the OG (needed for PREQUAL O(L^2) memory).
    ch_spptree_fas = spptree_og_map.merge(ORTHOFINDER_MCL_ALL.out.spptree_fas.flatten())
        .map { meta, fasta ->
            def text = fasta.text
            def n_seq = text.count('>')
            def cur = 0; def maxL = 0
            text.eachLine { line ->
                if (line.startsWith('>')) { if (cur > maxL) maxL = cur; cur = 0 }
                else { cur += line.trim().length() }
            }
            if (cur > maxL) maxL = cur
            [meta + [n_seq: n_seq, max_len: maxL], fasta]
        }
    ch_genetree_fas = genetree_og_map.merge(ORTHOFINDER_MCL_ALL.out.genetree_fas.flatten())
        .map { meta, fasta ->
            def text = fasta.text
            def n_seq = text.count('>')
            def cur = 0; def maxL = 0
            text.eachLine { line ->
                if (line.startsWith('>')) { if (cur > maxL) maxL = cur; cur = 0 }
                else { cur += line.trim().length() }
            }
            if (cur > maxL) maxL = cur
            [meta + [n_seq: n_seq, max_len: maxL], fasta]
        }

    // --test mode: keep only gene families that contain ALL species in the dataset.
    // This dramatically reduces the number of OGs for quick end-to-end smoke tests.
    if (params.test_run) {
        ch_total_species = renamed_prots.collect { it[0].id }
            .map { ids -> ids.unique().size() }

        ch_spptree_fas = ch_spptree_fas
            .combine(ch_total_species)
            .filter { meta, fasta, n_spp ->
                def species = fasta.text.readLines()
                    .findAll { it.startsWith('>') }
                    .collect { it.substring(1).replaceFirst(/_[^_]+$/, '') }
                    .toSet()
                species.size() >= n_spp
            }
            .map { meta, fasta, n_spp -> [meta, fasta] }

        ch_genetree_fas = ch_genetree_fas
            .combine(ch_total_species)
            .filter { meta, fasta, n_spp ->
                def species = fasta.text.readLines()
                    .findAll { it.startsWith('>') }
                    .collect { it.substring(1).replaceFirst(/_[^_]+$/, '') }
                    .toSet()
                species.size() >= n_spp
            }
            .map { meta, fasta, n_spp -> [meta, fasta] }
    }

    emit:
    spptree_fas   = ch_spptree_fas
    genetree_fas  = ch_genetree_fas
    inflation_dir = ORTHOFINDER_MCL_ALL.out.inflation_dir
}
