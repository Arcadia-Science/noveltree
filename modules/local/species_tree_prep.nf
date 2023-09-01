process SPECIES_TREE_PREP {
    tag "Preparing for species tree inference."
    label 'process_high'

    container "${ workflow.containerEngine == 'docker' ?
        'ubuntu:20.04': '' }"

    publishDir(
        path: "${params.outdir}/species_tree_prep",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    file genetrees  // Output from IQ-tree: filepaths to gene family trees with metadata
    file alignments // Output from ClipKit: filepaths to trimmed alignments with metadata
    val family_set  // String indicating whether these are gene families intended for SpeciesRax or GeneRax

    output:
    path "*gene_family_trees.txt" , emit: treefile
    path "*speciesrax_map.link"   , emit: speciesrax_map, optional: true
    path "*generax_map.link"      , emit: generax_map, optional: true
    path "asteroid_map.link"      , emit: asteroid_map, optional: true
    path "*orthogroup.families"   , emit: families
    path "*.family"               , emit: per_gene_family, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Take as input the trimmed "core" orthogroup MSAs from clipkit and
    # corresponding inferred gene family trees from IQ-tree...
    # Output a number of files necessary for downstream species tree
    # inference with both Asteroid and SpeciesRax, as well as gene-tree
    # species-tree reconciliation and gene-family duplication/transfer/loss
    # rates with GeneRax.

    # From these inputs we need to make:
    # 1) A file that contains each gene-family tree in newick
    #    format, one per line
    # 2) A tab separated species-gene mapping file, one gene-species pair per line
    #    - Two versions need to be made:
    #       a) Asteroid: a single file, with all mappings
    #       b) GeneRax: one file per gene-family
    # 3) A "families" file that corresponds each gene-family tree, multiple
    #    sequence alignment, and species-gene map.

    # Begin by creating the gene-tree file. Very simple.
    cat ./*treefile > ${family_set}_gene_family_trees.txt

    # Now, create the Species/GeneRax mapping files, and concatenate for Asteroid.
    # We will concurrently populate the families file.
    echo "[FAMILIES]" > ${family_set}_orthogroup.families
    for msa in \$(ls ./*fa)
    do
        # Get the OG name
        og=\$(echo \$msa | cut -f1 -d"_")
        tree=\$(ls \${og}*.treefile)

        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        grep ">" \${og}* | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > \${og}_${family_set}_map.link
        rm prot && rm spp

        if [[ ${family_set} == "generax" ]]; then
            # We need one family file per gene family for the per-family generax analysis.
            # Make this.
            echo "[FAMILIES]" > \${og}.family
            echo "- \${og}" >> \${og}.family
            echo "starting_gene_tree = \${tree}" >> \${og}.family
            echo "mapping = \${og}_${family_set}_map.link" >> \${og}.family
            echo "alignment = \$msa" >> \${og}.family
            echo "subst_model = LG+G4+F" >> \${og}.family
            sed -i 's|\\./||g' \${og}.family

            # Populate the families file for this gene family for the
            # Per-species analysis with generax
            # We will be using LG+G4+F for all gene families
            echo "- \${og}" >> ${family_set}_orthogroup.families
            echo "starting_gene_tree = \${og}_reconciled_gft.newick" >> ${family_set}_orthogroup.families
            echo "mapping = \${og}_${family_set}_map.link" >> ${family_set}_orthogroup.families
            echo "alignment = \$msa" >> ${family_set}_orthogroup.families
            echo "subst_model = LG+G4+F" >> ${family_set}_orthogroup.families
        else
            # Populate the families file for this gene family for the
            # analysis with SpeciesRax
            # We will be using LG+G4+F for all gene families
            echo "- \${og}" >> ${family_set}_orthogroup.families
            echo "starting_gene_tree = \${tree}" >> ${family_set}_orthogroup.families
            echo "mapping = \${og}_${family_set}_map.link" >> ${family_set}_orthogroup.families
            echo "alignment = \$msa" >> ${family_set}_orthogroup.families
            echo "subst_model = LG+G4+F" >> ${family_set}_orthogroup.families
        fi
    done

    # clean up the families file a bit
    sed -i 's|\\./||g' ${family_set}_orthogroup.families

    # Now concatenate the GeneRax maps for input to Asteroid
    if [[ ${family_set} == "generax" ]]
    then
      cat *${family_set}_map.link > asteroid_map.link
    fi
    """
}
