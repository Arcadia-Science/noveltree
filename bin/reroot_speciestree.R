#!/usr/bin/env Rscript
# Load necessary libraries
library(ape)
library(phytools)

reroot_tree <- function(asteroid_tree, outgroups){
    # Read in the newick tree file
    tree <- ape::unroot(ape::read.tree(asteroid_tree))
    
    # Reroot the tree
    if(length(outgroups) > 1){
        tree$node.label <- 1:tree$Nnode
        mrca <- ape::getMRCA(tree, outgroups)
        tree <- phytools::reroot(tree, node = mrca, position = 0.5)
    } else {
        tree$node.label <- 1:tree$Nnode
        parent_node <- 
            tree$edge[,1][which(tree$edge[,2] == match(outgroups, tree$tip.label))]
        tree <- phytools::reroot(tree, node = parent_node, position = 0.5)
    }
    
    # Write the rerooted tree to a file
    ape::write.tree(tree, file = "asteroid_rooted.bestTree.newick")
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    cat("Usage: Rscript reroot_speciestree.R <unrooted_speciestree> <outgroups>\n")
} else {
    # Read in the unrooted speciestree inferred with asteroid and the string 
    # containing the name(s) of outgroup species used to manually root it from
    # command line arguments
    asteroid_tree <- args[1]
    outgroups <- args[2]
    outgroups <- unlist(strsplit(outgroups, ","))

    # Call the function with the input and output file names
    reroot_tree(asteroid_tree, outgroups)
}
