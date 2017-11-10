# BIN-PHYLO-Tree_analysis.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PHYLO-Tree_analysis unit.
#
# Version:  1.0.1
#
# Date:     2017  10  31
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0.1  Wrong section heading
#           1.0    First 2017 version
#           0.1    First code copied from 2016 material.
#
#
# TODO:
#
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                       Line
#TOC> -------------------------------------------
#TOC>   1        ___Section___                 38
#TOC>   2        Tree Analysis                 77
#TOC>   2.1      Rooting Trees                136
#TOC>   2.2      Rotating Clades              182
#TOC>   2.3      Computing tree distances     229
#TOC>
#TOC> ==========================================================================


# =    1 Preparation and Tree Plot  ============================================


if (!require(Rphylip, quietly=TRUE)) {
  install.packages("Rphylip")
  library(Rphylip)
}
# Package information:
#  library(help = Rphylip)       # basic information
#  browseVignettes("Rphylip")    # available vignettes
#  data(package = "Rphylip")     # available datasets



# Read the species tree that you have created at the phyloT Website:
fungiTree <- read.tree("fungiTree.txt")

plot(fungiTree)

# The tree produced by phyloT contains full length species names, but it would
# be more convenient if it had bicodes instead.
str(fungiTree)

# The species names are in a vector $tip.label of this list. We can use bicode()
# to shorten them - but note that they have underscores as word separators. Thus
# we will use gsub("-", " ", ...) to replace the underscores with spaces.

for (i in seq_along(fungiTree$tip.label)) {
  fungiTree$tip.label[i] <- biCode(gsub("_", " ", fungiTree$tip.label[i]))
}

# Plot the tree
plot(fungiTree, cex=1.0, root.edge=TRUE, no.margin=TRUE)
nodelabels(text=orgTree$node.label, cex=0.6, adj=0.2, bg="#D4F2DA")
# Note that you can use the arrow buttons in the menu above the plot to scroll
# back to plots you have created earlier - so you can reference back to the
# species tree.


# =    2  Tree Analysis  =======================================================


# 1.1  Visualizing your tree
# The trees that are produced by Rphylip are stored as an object of class
# "phylo". This is a class for phylogenetic trees that is widely used in the
# community, practically all R phylogenetics packages will options to read and
# manipulate such trees. Outside of R, a popular interchange format is the
# Newick_format that you have seen above. It's easy to output your calculated
# trees in Newick format and visualize them elsewhere.

# The "phylo" class object is one of R's "S3" objects and methods to plot and
# print it have been defined with the Rphylip package, and the package ape that
# Rphylip has loaded. You can simply call plot(<your-tree>) and R knows what to
# do with <your-tree> and how to plot it. The underlying function is
# plot.phylo(), and documentation for its many options can by found by typing:

?plot.phylo

# We load the APSES sequence tree that you produced in the
# BIN-PHYLO-Tree_building unit:
load(file = "APSEStreeRproml.RData")

plot(apsTree) # default type is "phylogram"
plot(apsTree, type="unrooted")
plot(apsTree, type="fan", no.margin = TRUE)

# rescale to show all of the labels:
# record the current plot parameters by assigning them to a variable ...
(tmp <- plot(apsTree, type="fan", no.margin = TRUE, plot=FALSE))
# ... and adjust the plot limits for a new plot:
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

# ==   2.1  Rooting Trees  =====================================================

# In order to analyse the tree, it is helpful to root it first and reorder its
# clades. Contrary to documentation, Rproml() returns an unrooted tree.

is.rooted(apsTree)

# You can root the tree with the command root() from the "ape" package. ape is
# automatically installed and loaded with Rphylip.

plot(apsTree)

# add labels for internal nodes and tips
nodelabels(cex=0.5, frame="circle")
tiplabels(cex=0.5, frame="rect")

# The outgroup of the tree is tip "11" in my sample tree, it may be a different
# number in yours. Substitute the correct node number below for "outgroup".
apsTree <- root(apsTree, outgroup = 11, resolve.root = TRUE)
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


# ==   2.2  Rotating Clades  ===================================================

# To interpret the tree, it is useful to rotate the clades so that they appear
# in the order expected from the cladogram of species.

# We can either rotate around individual internal nodes ...
layout(matrix(1:2, 1, 2))
plot(apsTree, no.margin=TRUE, root.edge=TRUE)
nodelabels(node=17, cex=0.7, bg="#ff8866")
plot(rotate(apsTree, node=17), no.margin=TRUE, root.edge=TRUE)
nodelabels(node=17, cex=0.7, bg="#88ff66")
# Note that the species at the bottom of the clade descending from node
# 17 is now plotted at the top.
layout(matrix(1), widths=1.0, heights=1.0)

# ... or we can plot the tree so it corresponds as well as possible to a
# predefined tip ordering. Here we use the ordering that phyloT has returned
# for the species tree.

# (Nb. we need to reverse the ordering for the plot. This is why we use the
# expression [nOrg:1] below instead of using the vector directly.)

nOrg <- length(apsTree$tip.label)

layout(matrix(1:2, 1, 2))
plot(fungiTree,
     no.margin=TRUE, root.edge=TRUE)
nodelabels(text=fungiTree$node.label, cex=0.5, adj=0.2, bg="#D4F2DA")

plot(rotateConstr(apsTree, apsTree$tip.label[nOrg:1]),
     no.margin=TRUE, root.edge=TRUE)
add.scale.bar(length=0.5)
layout(matrix(1), widths=1.0, heights=1.0)

# Task: Study the two trees and consider their similarities and differences.
#         What do you expect? What do you find? Note that this is not a "mixed"
#         gene tree yet, since it contains only a single gene for the species
#         we considered. All of the branch points in this tree are speciation
#         events. Thus the gene tree should have the same topology as the
#         species tree. Does it? Are the differences important? How many
#         branches would you need to remove and reinsert elsewhere to get the
#         same topology as the species tree?

# In order to quantiofy how different these tow trees are, we need to compute
# tree distances.


# ==   2.3  Computing tree distances  ==========================================


# Many superb phylogeny tools are contributed by the phangorn package.

if (!require(phangorn, quietly=TRUE)) {
  install.packages("phangorn")
  library(phangorn)
}
# Package information:
#  library(help = phangorn)       # basic information
#  browseVignettes("phangorn")    # available vignettes
#  data(package = "phangorn")     # available datasets

# To compare two trees, they must have the same tip labels. We delete "MBP1_" or
# "KILA_" from the existing tip labels in a copy of our APSES domain tree.
apsTree2 <- apsTree
apsTree2$tip.label <- gsub("(MBP1_)|(KILA_)", "", apsTree2$tip.label)

# phangorn provides several functions to compute tree-differences (and there
# is a _whole_ lot of theory on how to compare trees). treedist() returns the
# "symmetric difference"
treedist(fungiTree, apsTree2, check.labels = TRUE)

# Numbers. What do they mean? How much more similar is our apsTree to the
# (presumably) ground truth of fungiTree than a random tree would be?
# The ape package (which was loaded with RPhylip) provides the function rtree()
# to compute random trees.

rtree(n = length(apsTree2$tip.label),  # number of tips
      rooted = TRUE,                   # we rooted the tree above,
                                       #  and fungiTree is rooted anyway
      tip.label = apsTree2$tip.label,  # use the apsTree2 labels
      br = NULL)                       # don't generate branch lengths since
                                       #   fungiTree has none, so we can't
                                       #   compare them anyway.

# Let's compute some random trees this way, calculate the distances to
# fungiTree, and then compare the values we get for apsTree2:

set.seed(112358)
N <- 10000  # takes about 15 seconds
myTreeDistances <- matrix(numeric(N * 2), ncol = 2)
colnames(myTreeDistances) <- c("symm", "path")

for (i in 1:N) {
  xTree <- rtree(n = length(apsTree2$tip.label),
                 rooted = TRUE,
                 tip.label = apsTree2$tip.label,
                 br = NULL)
  myTreeDistances[i, ] <- treedist(fungiTree, xTree)
}

table(myTreeDistances[, "symm"])

(symmObs <- treedist(fungiTree, apsTree2)[1])

# Random events less-or-equal to observation, divided by total number of
# events gives us the empirical p-value.
cat(sprintf("\nEmpirical p-value for symmetric diff. of observed tree is %1.4f\n",
            (sum(myTreeDistances[ , "symm"] <= symmObs) + 1) / (N + 1)))

hist(myTreeDistances[, "path"],
     col = "aliceblue",
     main = "Distances of random Trees to fungiTree")
(pathObs <- treedist(fungiTree, apsTree2)[2])
abline(v = pathObs, col = "chartreuse")

# Random events less-or-equal to observation, divided by total number of
# events gives us the empirical p-value.
cat(sprintf("\nEmpirical p-value for path diff. of observed tree is %1.4f\n",
            (sum(myTreeDistances[ , "path"] <= symmObs) + 1) / (N + 1)))

# Indeed, our apsTree is _very_ much more similar to the species tree than
# we would expect by random chance.

# What do we gain from that analysis? Analyzing the tree we get from a single
# gene of orthologous sequences is a positive control in our computational
# experiment. If these genes are indeed orthologues, a correct tree-building
# program ought to give us a tree that exactly matches the species tree.
# Evaluating how far off we are from the known correct result gives us a way to
# validate our workflow and our algorithm. If we can't get that right, we can't
# expect to get "real" data right either. Employing such positive controls in
# every computational experiment is essential for research. Not doing so is
# Cargo Cult Bioinformatics.


# [END]
