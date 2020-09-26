# tocID <- "BIN-PHYLO-Tree_analysis.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PHYLO-Tree_analysis unit.
#
# Version:  1.2
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2020 updates. Deprecate iTol and use taxize:: instead.
#                  Rewrite of tip re-ordering. Better handling of
#                  messages. pBar() for randomization.
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.0.2  Typo in variable name, style changes
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
#TOC>   Section  Title                              Line
#TOC> --------------------------------------------------
#TOC>   1        Preparation and Tree Plot            50
#TOC>   2        SPECIES REFERENCE TREE               66
#TOC>   3        Tree Analysis                       117
#TOC>   3.1        Rooting Trees                     177
#TOC>   3.2        Rotating Clades                   222
#TOC>   3.3        Computing tree distances          309
#TOC> 
#TOC> ==========================================================================


# =    1  Preparation and Tree Plot  ===========================================


if (! requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
# Package information:
#  library(help = ape)       # basic information
#  browseVignettes("ape")    # available vignettes
#  data(package = "ape")     # available datasets

# We change the graphics parameters from time to time, let's define the
# default so we can recreate a sane state:
dev.off()
PAR <- par()

# =    2  SPECIES REFERENCE TREE  ==============================================

# Before we do any kind of phylogenetic analysis of genes from several species,
# we MUST have a reference tree of the taxonomic relationships in hand. This
# context is absolutely required for the interpretation of our tree.

# We have the tax-ids in our database, and the NCBI has the species tree - we just need some way to extract the subtree that corresponds to our taxons of interest. Here's how to use the taxize:: package.

if (! requireNamespace("taxize", quietly = TRUE)) {
  install.packages("taxize")
}
# Package information:
#  library(help   = taxize)       # basic information
#  browseVignettes("taxize")    # available vignettes
#  data(package  = "taxize")     # available datasets

( mySOI <- c(myDB$taxonomy$ID, "83333") )
myClass <- taxize::classification(mySOI, db = "ncbi")
str(myClass)

myClass[[1]]

fungiTree <- taxize::class2tree(myClass, check = TRUE)
plot(fungiTree)

# The tree produced by taxize:: contains full length species names,
# but it would be more convenient if it had bicodes instead. Also, the actual
# tree is only part of the list(), which will cause problems later:
str(fungiTree)

# we therefor simplify
fungiTree <- fungiTree$phylo
str(fungiTree)

# The species names are in a vector $phylo$tip.label of this list.
# We can use biCode() to shorten them.
fungiTree$tip.label <- biCode(fungiTree$tip.label)

# Plot the tree
nSP <- length(fungiTree$tip.label)
plot(fungiTree, cex = 0.8, root.edge = TRUE, no.margin = TRUE)
text(-1, nSP - 0.5, "Species Tree:\nFungi", pos = 4)
ape::nodelabels(text = fungiTree$node.label,
                cex = 0.6,
                adj = 0.2,
                bg = "#D4F2DA")
# Note that you can use the arrow buttons in the menu above the plot pane to
# scroll back to plots you have created earlier - so you can reference back to
# this species tree in your later analysis.


# =    3  Tree Analysis  =======================================================


# 1.1  Visualizing your tree
# The trees that are produced by Rphylip are stored as an object of class
# "phylo". This is a class for phylogenetic trees that is widely used in the
# community, practically all R phylogenetics packages will options to read and
# manipulate such trees. Outside of R, a popular interchange format is the
# Newick_format that you have seen above. It's easy to output your calculated
# trees in Newick format and visualize them elsewhere.

# The "phylo" class object is one of R's "S3" objects and methods to plot and
# print it have been defined with the Rphylip package, and in ape. You can
# simply call plot(<your-tree>) and R knows what to do with <your-tree> and how
# to plot it. The underlying function is plot.phylo(), and documentation for its
# many options can by found by typing:

?plot.phylo

# We load the APSES sequence tree that you produced in the
# BIN-PHYLO-Tree_building unit:
apsTree <- readRDS(file = "data/APSEStreeRproml.rds")

plot(apsTree) # default type is "phylogram"
plot(apsTree, type = "unrooted")
plot(apsTree, type = "fan", no.margin = TRUE)

# rescale to show all of the labels:
# record the current plot parameters by assigning them to a variable ...
(tmp <- plot(apsTree, type="fan", no.margin = TRUE, plot=FALSE))
# ... and adjust the plot limits for a new plot:
plot(apsTree,
     type = "fan",
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
ape::nodelabels()
ape::edgelabels()
ape::tiplabels()

# show the number of nodes, edges and tips
ape::Nnode(apsTree)
ape::Nedge(apsTree)
ape::Ntip(apsTree)

par(PAR)   # reset graphics state

# Finally, write the tree to console in Newick format
ape::write.tree(apsTree)

# ==   3.1  Rooting Trees  =====================================================

# In order to analyse the tree, it is helpful to root it first and reorder its
# clades. Contrary to documentation, Rproml() returns an unrooted tree.

ape::is.rooted(apsTree)

# You can root the tree with the command root() from the "ape" package.

plot(apsTree)

# add labels for internal nodes and tips
ape::nodelabels(cex = 0.5, frame = "circle")
ape::tiplabels(cex = 0.5, frame = "rect")

# The outgroup of the tree (KILA ESCCO) is tip "11" in my sample tree, it may be a different
# number in yours. Substitute the correct node number below for "outgroup".
apsTree <- ape::root(apsTree, outgroup = 11, resolve.root = TRUE)
plot(apsTree)
ape::is.rooted(apsTree)

# This tree _looks_ unchanged, beacuse when the root trifurcation was resolved,
# an edge of length zero was added to connect the MRCA (Most Recent Common
# Ancestor) of the ingroup.

# The edge lengths are stored in the phylo object:
apsTree$edge.length

# ... and you can assign a small arbitrary value to the edge
# to show how it connects to the tree without having an
# overlap.
apsTree$edge.length[1] <- 0.1
plot(apsTree, cex = 0.7)
ape::nodelabels(text = "MRCA", node = 12, cex = 0.5, adj = 0.1, bg = "#ff8866")


# This procedure does however not assign an actual length to a root edge, and
# therefore no root edge is visible on the plot. Why? , you might ask. I ask
# myself that too. We'll just add a length by hand.

apsTree$root.edge <- mean(apsTree$edge.length) * 1.5
plot(apsTree, cex = 0.7, root.edge = TRUE)
ape::nodelabels(text = "MRCA", node = 12, cex = 0.5, adj = 0.8, bg = "#ff8866")


# ==   3.2  Rotating Clades  ===================================================

# To interpret the tree, it is useful to rotate the clades so that they appear
# in the order expected from the cladogram of species.

# We can either rotate around individual internal nodes ...
layout(matrix(1:2, 1, 2))
plot(apsTree, no.margin = TRUE, root.edge = TRUE)
ape::nodelabels(node = 13, cex = 0.7, bg = "#ff8866")
plot(ape::rotate(apsTree, node = 13), no.margin = TRUE, root.edge = TRUE)
ape::nodelabels(node = 13, cex = 0.7, bg = "#88ff66")
# Note that the species at the bottom of the clade descending from node
# 17 is now plotted at the top.

par(PAR)   # reset graphics state

# ... or we can rearrange the tree so it corresponds as well as possible to a
# predefined tip ordering. Here we use the ordering that taxize:: has inferred
# from the NCBI taxonomic classification.

nOrg <- length(apsTree$tip.label)

plot(fungiTree,
     no.margin = FALSE, root.edge = TRUE)
ape::nodelabels(text = fungiTree$node.label,
                cex = 0.5,
                adj = 0.2,
                bg = "#D4F2DA")

# These are the fungi tree tips ...
fungiTree$tip.label
# ... and their order is determined by the edge-list that is stored in
fungiTree$edge
# which edges join the tips?
ape::tiplabels(cex = 0.5, frame = "rect")
# as you can see, the tips (range [1:nOrg] ) are in column 2 and they are
# ordered from bottom to top.
# And each tip number is the index of the species in the tip.label vector. So we can take column 2, subset it, and use it to get a list of species in the order of the tree ...

sel <- fungiTree$edge[ , 2 ] <= nOrg
( oSp <- fungiTree$tip.label[fungiTree$edge[sel , 2 ]] )

# Now, here are the genes of the apsTree tips ...
apsTree$tip.label

# ... and the "constraint"  we need for reordering, according to the help page
# of ape::rotateConstr(), is "a vector specifying the order of the tips as they
# should appear (from bottom to top)". Thus we need to add the "MBP1_" prefix to our vector
oSp <- gsub("^", "MBP1_", oSp)
( oSp <- gsub("MBP1_ESSCO", "KILA_ESCCO", oSp) )

# Then we can plot the two trees to compare: the fungi- tree
par(PAR)   # reset graphics state
layout(matrix(1:2, 1, 2))
plot(fungiTree,
    no.margin = TRUE,
     root.edge = TRUE)
ape::nodelabels(text = fungiTree$node.label,
                cex = 0.5,
                adj = 0.2,
                bg = "#D4F2DA")

# and the re-organized apsesTree ...
plot(ape::rotateConstr(apsTree, constraint = oSp[]),
     no.margin = TRUE,
     root.edge = TRUE)

par(PAR)   # reset graphics state

# As you can see, the reordering is not perfect, since the topologies are
# different, mostly due to the unresolved nodes in the reference tree. One
# could play with that ...


# Task: Study the two trees and consider their similarities and differences.
#         What do you expect? What do you find? Note that this is not a "mixed"
#         gene tree yet, since it contains only a single gene for the species
#         we considered. All of the branch points in this tree are speciation
#         events. Thus the gene tree should have the same topology as the
#         species tree. Does it? Are the differences important? How many
#         branches would you need to remove and reinsert elsewhere to get the
#         same topology as the species tree?

# In order to quantify how different these two trees are, we need to compute
# tree distances.


# ==   3.3  Computing tree distances  ==========================================


# Many superb phylogeny tools are contributed by the phangorn package.

if (! requireNamespace("phangorn", quietly = TRUE)) {
  install.packages("phangorn")
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
phangorn::treedist(fungiTree, apsTree2, check.labels = TRUE)

# Numbers. What do they mean? How much more similar is our apsTree to the
# (presumably) ground truth of fungiTree than a random tree would be?
# The ape package provides the function rtree()
# to compute random trees.

ape::rtree(n = length(apsTree2$tip.label), # number of tips
          rooted = TRUE,                   # we rooted the tree above,
                                           #  and fungiTree is rooted anyway
          tip.label = apsTree2$tip.label,  # use the apsTree2 labels
          br = NULL)                       # don't generate branch lengths since
                                           #   fungiTree has none, so we can't
                                           #   compare them anyway.

# (Note the warning message about non-binary trees; we'll suppress that later
#  by wrapping the function call in supressMessages(); we don't want to
#  print it 10,000 times :-)


# Let's compute some random trees this way, calculate the distances to
# fungiTree, and then compare the values we get for apsTree2. The random
# trees are provided by ape::rtree().

N <- 10000  # takes about 15 seconds, and we'll use the pBar function,
            # defined in .utilities.R  to keep track of where we are at:
myTreeDistances <- matrix(numeric(N * 2), ncol = 2)
colnames(myTreeDistances) <- c("symm", "path")

set.seed(112358)
for (i in 1:N) {
  pBar(i, N)
  xTree <- ape::rtree(n = length(apsTree2$tip.label),
                      rooted = TRUE,
                      tip.label = apsTree2$tip.label,
                      br = NULL)
  myTreeDistances[i, ] <- suppressMessages(phangorn::treedist(fungiTree, xTree))
}
set.seed(NULL)                      # reset the random number generator

table(myTreeDistances[, "symm"])

( symmObs <- phangorn::treedist(fungiTree, apsTree2)[1] )

# Random events less-or-equal to observation, divided by total number of
# events gives us the empirical p-value.
cat(sprintf("\nEmpirical p-value for symmetric diff. of observed tree is %1.4f\n",
            (sum(myTreeDistances[ , "symm"] <= symmObs) + 1) / (N + 1)))

par(PAR)   # reset graphics state
hist(myTreeDistances[, "path"],
     col = "aliceblue",
     main = "Distances of random Trees to fungiTree")
(pathObs <- phangorn::treedist(fungiTree, apsTree2)[2])
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
