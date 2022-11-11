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
#TOC>   1        Packages                                55
#TOC>   2        PhyML online server workflow            68
#TOC> 
#TOC> ==========================================================================


# ==============================================================================
#
#      I N   P R O G R E S S
#
#   I am removing PHYYLIP from the course material since it is no longer
#   maintained and RPhylip has recently been removed from CRAN. This script
#   is now functional although it will still be updated over the day.
#
# ==============================================================================


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

# This will install phangorn, as well as its dependency, the package "ape".

# =    2  PhyML online server workflow  ========================================

# read the multi-FASTA alignment
tmp <- phangorn::read.phyDat("data/APSESphyloSet.mfa",
                             format = "fasta",
                             type = "AA")

# Save ...
# phangorn::write.phyDat(tmp, file = "data/APSESphyloSet.phy")

#
# Submit the file to http://www.atgc-montpellier.fr/phyml/
#    1.  Navigate to the Website
#    2.  Upload "data/APSESphyloSet.phy"
#    3.  Use default parameters.
#    4.  Make sure to enter you eMail address to be notified of the results
#
# The computation may complete in a minute or so. It may also take longer.
#    5.  Download the .zip attachment to your results email and expand
#    6.  From the folder with results find "apsesphyloset_phy_phyml_tree.txt"
#          and move it to your data/ folder

# Confirm that you can read the tree and that it makes sense
apsTree <- ape::read.tree("data/apsesphyloset_phy_phyml_tree.txt")

plot(apsTree)

# If this did not work, ask for advice.



# [END]
