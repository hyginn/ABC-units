# tocID <- "BIN-PHYLO-Tree_building.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PHYLO-Tree_building unit.
#
# Version:  2.0
#
# Date:     2017 - 10  -  2022 - 11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.0    The PHYLIP era is over
#           1.2    deprecate save()/load() for saveRDS()/readRDS(); Mac:
#                  instructions to authorize proml.app
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#           1.0    First 2017 version
#           0.1    First code copied from 2016 material.
#
# Note:
#           This unit was originally developed as a workflow with the
#           groundbreaking maximum-likelihood methods that were developed
#           by Joe Felsenstein in Washington as part of the Phylip suite.
#           Unfortunately, Phylip is no longer actively maintained, and
#           although the code works as designed, its arkane user interface
#           has made it unsuitable for coursework by novices. For a while,
#           the Rphylip:: package provided a useable wrapper for
#           R-scripted analysis, but since that too has gone out of maintenance,
#           and the package was removed from CRAN in June 2022, I am, with
#           gratitude and respect, removing Phylip from the course. (2022-11)
#
#           cf. https://evolution.genetics.washington.edu/phylip.html
#
#
# TODO:
#           Add MrBayes
# https://cran.r-project.org/web/packages/phangorn/vignettes/IntertwiningTreesAndNetworks.html
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
#TOC>   Section  Title                                 Line
#TOC> -----------------------------------------------------
#TOC>   1        Packages                                60
#TOC>   2        PhyML online server workflow            78
#TOC>   2.1        .mfa to .phy                          83
#TOC>   2.2        Computing the tree                   100
#TOC>   2.3        Reading the tree back into R         114
#TOC>
#TOC> ==========================================================================


# =    1  Packages  ============================================================
#
#
if (! requireNamespace("phangorn", quietly = TRUE)) {
  install.packages("phangorn")
}
# Package information:
#  library(help = phangorn)       # basic information
#  browseVignettes("phangorn")    # available vignettes
#  data(package = "phangorn")     # available datasets

# This will install phangorn::, as well as its dependency, the package ape::.
# Here, we only use phangorn:: to read a multi-FASTA file and write a
# Phylip-formatted dataset that is suitable as input for tree-inference
# programs. But phangorn:: can do a lot more than that and even has its own
# maximum-likelihood tree-inference code: phangorn::pml().


# =    2  PhyML online server workflow  ========================================

# Workflow to create input that is suitable for an online version of PhyML.


# ==   2.1  .mfa to .phy  ======================================================

# In this workflow we reformat a multi-FASTA file into the .phy format that
# is used as input for tree inference by many different programs. The first
# line needs to specify the number of organisms and the number of "states"
# (sequence characters), subsequent lines contain the organism name and the
# data.

# Read the multi-FASTA alignment that we produced previously
tmp <- phangorn::read.phyDat("data/APSESphyloSet.mfa",
                             format = "fasta",
                             type = "AA")

# Write the alignment to disk in Phylip format
# phangorn::write.phyDat(tmp, file = "data/APSESphyloSet.phy")


# ==   2.2  Computing the tree  ================================================

# Submit the file to the Montpellier PhyML server
#    1.  Navigate to http://www.atgc-montpellier.fr/phyml/
#    2.  Upload "data/APSESphyloSet.phy"
#    3.  Use default parameters
#    4.  Make sure to enter your eMail address to be notified of the results
#
# The computation may complete in a minute or so. It may also take longer.
#    5.  Download the .zip attachment to your results email and expand
#    6.  From the folder with results find "apsesphyloset_phy_phyml_tree.txt"
#          and move it to your data/ folder


# ==   2.3  Reading the tree back into R  ======================================

# Confirm that you can read the tree and that it makes sense
apsTree <- ape::read.tree("data/apsesphyloset_phy_phyml_tree.txt")

plot(apsTree)

# If this did not work, ask for advice.



# [END]
