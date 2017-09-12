# BIN-ALI-Similarity.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-Similarity unit.
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

# = 1 Mutation Data matrix

# First, we install and load the Biostrings package.
if (!require(Biostrings, quietly=TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
  library(Biostrings)
}


# Biostrings contains mutation matrices and other useful datasets
data(package = "Biostrings")

# Let's load BLOSUM62
data(BLOSUM62)

# ... and see what it contains. (You've seen this before, right?)
BLOSUM62

# We can simply access values via the row/column names to look at the data
# for the questions I asked in the Assignment on the Wiki:
BLOSUM62["H", "H"]
BLOSUM62["S", "S"]

BLOSUM62["L", "K"]
BLOSUM62["L", "I"]


BLOSUM62["R", "W"]
BLOSUM62["W", "R"]   # the matrix is symmetric!




# = 1.1 <<<Subsection>>>



# = 1 Tasks




# [END]
