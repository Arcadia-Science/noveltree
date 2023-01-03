process SPECIES_TREE_PREP {
    tag "Preparing for species tree inference."
    label 'process_single'

    container "ubuntu:20.04"
    
    publishDir(
        path: "${params.outdir}/species_tree_prep",
        mode: params.publish_dir_mode,
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path genetrees   // Output from IQ-tree: filepaths to gene family trees with metadata
    path alignments // Output from ClipKit: filepaths to trimmed alignments with metadata

    output:
    path "gene-family-trees.txt" ,       emit: treefile
    path "*generax-map.link" ,           emit: generax_map
    path "asteroid-map.link" ,           emit: asteroid_map
    path "generax-orthogroup.families" , emit: families

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
    cat ./*treefile > gene-family-trees.txt

    # Now, create the GeneRax mapping files, and concatenate for Asteroid.
    # We will concurrently populate the families file.
    echo "[FAMILIES]" > generax-orthogroup.families
    for msa in \$(ls ./*fa)
    do
        # Get the OG name
        og=\$(echo \$msa | sed "s/-clipkit.fa//g")
        tree=\$(ls \${og}*.treefile)

        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species
        # and then the protein
        grep ">" \${og}* | sed "s/>//g"  | sed "s/.*://g" > prot
        sed "s/_[^_]*\$//" prot | sed "s/EP0*._//g" > spp
        paste prot spp > \${og}-generax-map.link
        rm prot && rm spp

        # Populate the families file for this gene family
        # We will be using LG+G4+F for all gene families
        echo "- \${og}" >> generax-orthogroup.families
        echo "starting_gene_tree = \${tree}" >> generax-orthogroup.families
        echo "mapping = \${og}-generax-map.link" >> generax-orthogroup.families
        echo "alignment = \$msa" >> generax-orthogroup.families
        echo "subst_model = LG+G4+F" >> generax-orthogroup.families
    done

    # clean up the families file a bit
    sed -i 's|\\./||g' generax-orthogroup.families

    # Now concatenate the maps for input to Asteroid
    cat *generax-map.link > asteroid-map.link
    """
}
