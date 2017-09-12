# RPR-Biostrings.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Biostrings unit.
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

# = 1 The Biostrings package

# First, we install and load the Biostrings package.
if (!require(Biostrings, quietly=TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
  library(Biostrings)
}

# This is a large collection of tools ...
help(package = "Biostrings")





# = 1.1 <<<Subsection>>>



# = 1 Tasks




# [END]
