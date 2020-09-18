# tocID <- "BIN-PHYLO-Data_preparation.R"
#
# ---------------------------------------------------------------------------- #
#  PATIENCE  ...                                                               #
#    Do not yet work wih this code. Updates in progress. Thank you.            #
#    boris.steipe@utoronto.ca                                                  #
# ---------------------------------------------------------------------------- #
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PHYLO-Data_preparation unit.
#
# Version:  1.1
#
# Date:     2017  10  -  2019  01
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
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
#TOC>   Section  Title                                     Line
#TOC> ---------------------------------------------------------
#TOC>   1        Preparations                                44
#TOC>   2        Fetching sequences                          76
#TOC>   3        Multiple Sequence Alignment                117
#TOC>   4        Reviewing and Editing Alignments           136
#TOC>   4.1        Masking workflow                         152
#TOC>
#TOC> ==========================================================================


# =    1  Preparations  ========================================================


# You need to reload your protein database, including changes that might have
# been made to the reference files. If you have worked with the prerequiste
# units, you should have a script named "makeProteinDB.R" that will create the
# myDB object with a protein and feature database. Ask for advice if not.
source("makeProteinDB.R")

# Load packages we need

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets


if (! requireNamespace("msa", quietly = TRUE)) {
  BiocManager::install("msa")
}
# Package information:
#  library(help = msa)       # basic information
#  browseVignettes("msa")  # available vignettes
#  data(package = "msa")   # available datasets


# =    2  Fetching sequences  ==================================================


# myDB contains the ten Mbp1 orthologues from the reference species and the Mbp1
# RBM for MYSPE. We will construct a phylogenetic tree from the proteins' APSES
# domains. You have annotated their ranges as a feature. The following code
# retrieves the sequences from myDB. You have seen similar code in other units.

sel <- grep("^MBP1_", myDB$protein$name)
(proNames <- myDB$protein$name[sel])
(proIDs <- myDB$protein$ID[sel])

(sel <- myDB$feature$ID[myDB$feature$name == "APSES fold"])
(fanIDs <- myDB$annotation$ID[myDB$annotation$proteinID %in% proIDs & # %in% !
                              myDB$annotation$featureID == sel])      #  ==  !
                                                                      # Why?
APSI <- character(length(fanIDs))

for (i in seq_along(fanIDs)) {
  sel   <- myDB$annotation$ID == fanIDs[i]  # get the feature row index
  proID <- myDB$annotation$proteinID[sel]   # get its protein ID
  start <- myDB$annotation$start[sel]       # get start ...
  end   <- myDB$annotation$end[sel]         # ... and end

  sel <- myDB$protein$ID == proID           # get the protein row index ...
                                            # ... and the sequence
  APSI[i] <- substring(myDB$protein$sequence[sel], start, end)
  names(APSI)[i] <- (myDB$protein$name[sel])
}

head(APSI)

# Let's add the E.coli Kila-N domain sequence as an outgroup, for rooting our
# phylogenetic tree (see the unit's Wiki page for details on the sequence).

APSI <- c(APSI,
"IDGEIIHLRAKDGYINATSMCRTAGKLLSDYTRLKTTQEFFDELSRDMGIPISELIQSFKGGRPENQGTWVHPDIAINLAQ")
names(APSI)[length(APSI)] <- "KILA_ESCCO"
tail(APSI)


# =    3  Multiple Sequence Alignment  =========================================

# This vector of sequences with named elements fulfills the requirements to be
# imported as a Biostrings object - an AAStringSet - which we need as input for
# the MSA algorithms in Biostrings.
#

APSESSet <- Biostrings::AAStringSet(APSI)
APSESMsa <- msa::msaMuscle(APSESSet, order = "aligned")

# Nb. msaMuscle() sometimes fails - reproducibly, but I am not sure why. If
# that happens in your case, just use msaClustalOmega() instead.

# inspect the alignment.
writeALN(APSESMsa)

# What do you think? Is this a good alignment for phylogenetic inference?


# =    4  Reviewing and Editing Alignments  ====================================


# Head back to the Wiki page for this unit and read up on the background
# first.

# Let's mask out all columns that have observations for
# less than 1/3 of the sequences in the dataset. This
# means they have more than round(nrow(msaSet) * (2/3))
# hyphens in a column.
#
# We take all sequences, split them into single
# characters, and put them into a matrix. Then we
# go through the matrix, column by column and decide
# whether we want to include that column.

# ==   4.1  Masking workflow  ==================================================

# get the length of the alignment
(lenAli <- APSESMsa@unmasked@ranges@width[1])

# initialize a matrix that can hold all characters
# individually
msaMatrix <- matrix(character(nrow(APSESMsa) * lenAli),
                    ncol = lenAli)

# assign the correct rownames
rownames(msaMatrix) <- APSESMsa@unmasked@ranges@NAMES
for (i in 1:nrow(APSESMsa)) {
  msaMatrix[i, ] <- unlist(strsplit(as.character(APSESMsa@unmasked[i]), ""))
}

# inspect the result
msaMatrix[1:7, 1:14]

# Now let's make a logical vector with an element for each column that selects
# which columns should be masked out.

# The number of hyphens in a column is easy to count. Consider:

    msaMatrix[ , 20]
    msaMatrix[ , 20] == "-"
sum(msaMatrix[ , 20] == "-")

# Thus filling our logical vector is simple:

# initialize a mask
colMask <- logical(ncol(msaMatrix))

# define the threshold for rejecting a column
limit <- round(nrow(APSESMsa) * (2/3))

# iterate over all columns, and write TRUE if there are less-or-equal to "limit"
# hyphens, FALSE if there are more - i.e. TRUE columns will be used fr analysis
# and FALSE columns will be rejected.
for (i in 1:ncol(msaMatrix)) {
  count <- sum(msaMatrix[ , i] == "-")
  colMask[i] <- count <= limit # TRUE if less-or-equal to limit, FALSE if not
}

# Inspect the mask
colMask

# How many positions are being kept?
sum(colMask)

cat(sprintf("We are masking %4.2f %% of alignment columns.\n",
            100 * (1 - (sum(colMask) / length(colMask)))))


# Next, we use colMask to remove the masked columns from the matrix
# in one step:
maskedMatrix <- msaMatrix[ , colMask]

# check:
ncol(maskedMatrix)

# ... then collapse each row of single characters back into a string ...
APSESphyloSet <- character()
for (i in 1:nrow(maskedMatrix)) {
  APSESphyloSet[i] <- paste(maskedMatrix[i, ], collapse="")
}
names(APSESphyloSet) <- rownames(maskedMatrix)

# inspect ...
writeALN(APSESphyloSet)

# As you see, we have removed a three residue insertion from MBP1_NEUCR, and
# several indels from the KILA_ESCCO outgroup sequence.


# We save the aligned, masked domains to a file in multi-FASTA format.
writeMFA(APSESphyloSet, myCon = "APSESphyloSet.mfa")




# [END]
