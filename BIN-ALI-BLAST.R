# BIN-ALI-BLAST.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-BLAST unit.
#
# Version:  1.0
#
# Date:     2017  10  23
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    First live version 2017.
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
#TOC>   Section  Title                         Line
#TOC> ---------------------------------------------
#TOC>   1        Preparations                    41
#TOC>   2        Defining the APSES domain       54
#TOC>   3        Executing the BLAST search      76
#TOC>   4        Analysing results               98
#TOC> 
#TOC> ==========================================================================


# =    1  Preparations  ========================================================

if (!require(Biostrings, quietly=TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
  library(Biostrings)
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets


# =    2  Defining the APSES domain  ===========================================

# Load your protein database
source("makeProteinDB.R")

# Get the APSES domain sequence for MBP1_MYSPE feature annotation. (You have
# entered this data in the BIN-ALI-Optimal_sequence_alignment unit.)

(proID <- myDB$protein$ID[myDB$protein$name == "MBP1_<MYSSPE>"]) # <<< EDIT
(ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"])
(fanID <- myDB$annotation$ID[myDB$annotation$proteinID == proID &
                               myDB$annotation$featureID == ftrID])
(start <- myDB$annotation$start[myDB$annotation$ID == fanID])
(end   <- myDB$annotation$end[myDB$annotation$ID == fanID])
(apses <- substr(myDB$protein$sequence[myDB$protein$ID == proID],
                 start,
                 end))

# The MYSPE "apses" sequence is the sequence that we will use for our reverse
# BLAST search.


# =    3  Executing the BLAST search  ==========================================

# The ./scripts/BLAST.R code defines two functions to access the BLAST interface
# through its Web API, and to parse results. Have a look at the script, then
# source it:

source("./scripts/BLAST.R")

# Use BLAST() to find the best match to the MYSPE APSES domain in Saccharomyces
# cerevisiae:

BLASThits <- BLAST(apses,                       # MYSPE APSES domain sequence
                  db = "refseq_protein",        # database to search in
                  nHits = 10,                   #
                  E = 0.01,                     #
                  limits = "txid559292[ORGN]")  # S. cerevisiae S288c


length(BLASThits)  # There should be at least one hit there. Ask for advice
                   # in case this step fails.


# =    4  Analysing results  ===================================================

# The BLAST.R script has defined a convenience function to parse BLAST
# alignments.

(topHit <- parseBLASTalignment(BLASThits, idx = 1))   # Parse the top hit

# What is the refseq ID of the top hit
topHit$accession

# If this is "NP_010227.1" you have confirmed the RBM of the MYSPE apses
# domain. If it is not, ask me for advice.





# [END]
