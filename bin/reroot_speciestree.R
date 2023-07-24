#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Extract arguments
asteroid_tree <- args[1]
outgroups <- args[2]
outgroups <- unlist(strsplit(outgroups, ","))

# Load necessary libraries
library(ape)
library(phytools)

# Read in the newick tree file
tree <- phytools::unroot(ape::read.tree(asteroid_tree))

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
