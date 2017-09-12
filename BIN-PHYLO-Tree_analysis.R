# BIN-PHYLO-Tree_analysis.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PHYLO-Tree_analysis unit.
#
# Version:  0.1
#
# Date:     2017  08  28
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.1    First code copied from 2016 material.

#
# TODO:
#
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...

# ==============================================================================

# = 1 ___Section___

# ==============================================================================
#        PART FIVE: Tree analysis
# ==============================================================================

# A Entrez restriction command
cat(paste(paste(c(myDB$taxonomy$ID, "83333"), "[taxid]", sep=""), collapse=" OR "))

# The Common Tree from NCBI
# Download the EDITED phyliptree.phy
commonTree <- read.tree("phyliptree.phy")

# Plot the tree
plot(commonTree, cex=1.0, root.edge=TRUE, no.margin=TRUE)
nodelabels(text=commonTree$node.label, cex=0.6, adj=0.2, bg="#D4F2DA")


# === Visualizing your tree ====================================================

# The trees that are produced by Rphylip are stored as an object of class
# "phylo". This is a class for phylogenetic trees that is widely used in the
# community, practically all R phylogenetics packages will options to read and
# manipulate such trees. Outside of R, a popular interchange format is the
# Newick_format that you have seen above. It's easy to output your calculated
# trees in Newick format and visualize them elsewhere.

# The "phylo" class object is one of R's "S3" objects and methods to plot and
# print it have been added to the system. You can simply call plot(<your-tree>)
# and R knows what to do with <your-tree> and how to plot it. The underlying
# function is plot.phylo(), and documentation for its many options can by found
# by typing:

?plot.phylo

plot(apsTree) # default type is "phylogram"
plot(apsTree, type="unrooted")
plot(apsTree, type="fan", no.margin = TRUE)

# rescale to show all of the labels:
# record the current plot parameters ...
tmp <- plot(apsTree, type="fan", no.margin = TRUE, plot=FALSE)
# ... and adjust the plot limits for a new plot
plot(apsTree,
     type="fan",
     x.lim = tmp$x.lim * 1.8,
     y.lim = tmp$y.lim * 1.8,
     cex = 0.8,
     no.margin = TRUE)

# Inspect the tree object
str(apsTree)
apsTree$tip.label
apsTree$edge
apsTree$edge.length

# show the node / edge and tip labels on a plot
plot(apsTree)
nodelabels()
edgelabels()
tiplabels()

# show the number of nodes, edges and tips
Nnode(apsTree)
Nedge(apsTree)
Ntip(apsTree)


# Finally, write the tree to console in Newick format
write.tree(apsTree)

# === Rooting Trees ============================================================

# In order to analyse the tree, it is helpful to root it first and reorder its
# clades. Contrary to documentation, Rproml() returns an unrooted tree.

is.rooted(apsTree)

# You can root the tree with the command root() from the "ape" package. ape is
# automatically installed and loaded with Rphylip.

plot(apsTree)

# add labels for internal nodes and tips
nodelabels(cex=0.5, frame="circle")
tiplabels(cex=0.5, frame="rect")

# The outgroup of the tree is tip "8" in my sample tree, it may be a different
# number in yours. Substitute the correct node number below for "outgroup".
apsTree <- root(apsTree, outgroup = 8, resolve.root = TRUE)
plot(apsTree)
is.rooted(apsTree)

# This tree _looks_ unchanged, beacuse when the root trifurcation was resolved,
# an edge of length zero was added to connect the MRCA (Most Recent Common
# Ancestor) of the ingroup.

# The edge lengths are stored in the phylo object:
apsTree$edge.length

# ... and you can assign a small arbitrary value to the edge
# to show how it connects to the tree without having an
# overlap.
apsTree$edge.length[1] <- 0.1
plot(apsTree, cex=0.7)
nodelabels(text="MRCA", node=12, cex=0.5, adj=0.1, bg="#ff8866")


# This procedure does however not assign an actual length to a root edge, and
# therefore no root edge is visible on the plot. Why? , you might ask. I ask
# myself that too. We'll just add a length by hand.

apsTree$root.edge <- mean(apsTree$edge.length) * 1.5
plot(apsTree, cex=0.7, root.edge=TRUE)
nodelabels(text="MRCA", node=12, cex=0.5, adj=0.8, bg="#ff8866")


# === Rotating Clades ==========================================================

# To interpret the tree, it is useful to rotate the clades so that they appear
# in the order expected from the cladogram of species.

# We can either rotate around individual internal nodes:
layout(matrix(1:2, 1, 2))
plot(apsTree, no.margin=TRUE, root.edge=TRUE)
nodelabels(node=17, cex=0.7, bg="#ff8866")
plot(rotate(apsTree, node=17), no.margin=TRUE, root.edge=TRUE)
nodelabels(node=17, cex=0.7, bg="#88ff66")
layout(matrix(1), widths=1.0, heights=1.0)

# ... or we can plot the tree so it corresponds as well as possible to a
# predefined tip ordering. Here we use the ordering that NCBI Global Tree
# returns for the reference species - we have used it above to make the vector
# apsMbp1Names. You inserted your YFO name into that vector - but you should
# move it to its correct position in the cladogram.

# (Nb. we need to reverse the ordering for the plot. This is why we use the
# expression [nOrg:1] below instead of using the vector directly.)

nOrg <- length(apsTree$tip.label)

layout(matrix(1:2, 1, 2))
plot(commonTree,
     no.margin=TRUE, root.edge=TRUE)
nodelabels(text=commonTree$node.label, cex=0.5, adj=0.2, bg="#D4F2DA")

plot(rotateConstr(apsTree, apsTree$tip.label[nOrg:1]),
     no.margin=TRUE, root.edge=TRUE)
add.scale.bar(length=0.5)
layout(matrix(1), widths=1.0, heights=1.0)

# Study the two trees and consider their similarities and differences. What do
# you expect? What do you find?
#

# Print the two trees on one sheet of paper, write your name and student number,
# and bring it to class as your deliverable for this assignment. Also write two
# or three sentences about if/how the gene tree matches the species tree or not.


# = 1 Tasks




# [END]
