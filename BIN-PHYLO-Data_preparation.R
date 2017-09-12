# BIN-PHYLO-Data_preparation.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PHYLO-Data_preparation unit.
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
#        PART ONE: Choosing sequences
# ==============================================================================

# Start by loading libraries. You already have the packages installed.
library(Biostrings)
library(msa)
library(stringr)

# What is the latest version of myDB that you have saved?
list.files(pattern = "myDB.*")

# ... load it (probably myDB.05.RData - if not, change the code below).
load("myDB.05.RData")

# The database contains the ten Mbp1 orthologues from the reference species
# and the Mbp1 RBM for YFO.
#
# We will construct a phylogenetic tree from the proteins' APSES domains.
# You have annotated their ranges as a feature.

# Collect APSES domain sequences from your database. The function
# dbGetFeatureSequence() retrieves the sequence that is annotated for a feature
# from its start and end coordinates. Try:

dbGetFeatureSequence(myDB, "MBP1_SACCE", "APSES fold")

# Lets put all APSES sequences into a vector:
APSESnames <- myDB$protein$name[grep("^MBP1_", myDB$protein$name)]
APSES <- character(length(APSESnames))

for (i in 1:length(APSESnames)) {
  APSES[i] <- dbGetFeatureSequence(myDB, APSESnames[i], "APSES fold")
}

# Let's name the rows of our vector with the BiCode part of the protein name.
# This is important so we can keep track of which sequence is which. We use the
# gsub() funcion to substitute "" for "MBP1_", thereby deleting this prefix.
names(APSES) <- gsub("^MBP1_", "", APSESnames)

# inspect the result: what do you expect? Is this what you expect?
head(APSES)

# Let's add the E.coli Kila-N domain sequence as an outgroup, for rooting our
# phylogegetic tree (see the Assignment Course Wiki page for details on the
# sequence).

APSES[length(APSES) + 1] <-
  "IDGEIIHLRAKDGYINATSMCRTAGKLLSDYTRLKTTQEFFDELSRDMGIPISELIQSFKGGRPENQGTWVHPDIAINLAQ"
names(APSES)[length(APSES)] <- "ESCCO"


# ==============================================================================
#        PART TWO: Multiple sequence alignment
# ==============================================================================

# This vector of sequences with named elements fulfills the requirements to be
# imported as a Biostrings object - an AAStringSet - which we need as input for
# the MSA algorithms in Biostrings.
#

APSESSeqSet <- AAStringSet(APSES)

APSESMsaSet <- msaMuscle(APSESSeqSet, order = "aligned")

# inspect the alignment.
writeSeqSet(APSESMsaSet, format = "ali")


# What do you think? Is this a good alignment for phylogenetic inference?

# ==============================================================================
#        PART THREE: reviewing and editing alignments
# ==============================================================================

# Head back to the assignment 7 course wiki page and read up on the background
# first.
#



# Let's mask out all columns that have observations for
# less than 1/3 of the sequences in the dataset. This
# means they have more than round(nrow(msaSet) * (2/3))
# hyphens in a column.
#
# We take all sequences, split them into single
# characters, and put them into a matrix. Then we
# go through the matrix, column by column and decide
# whether we want to include that column.

# Step 1. Go through this by hand...

# get the length of the alignment
lenAli <- APSESMsaSet@unmasked@ranges@width[1]

# initialize a matrix that can hold all characters
# individually
msaMatrix <- matrix(character(nrow(APSESMsaSet) * lenAli),
                    ncol = lenAli)

# assign the correct rownames
rownames(msaMatrix) <- APSESMsaSet@unmasked@ranges@NAMES
for (i in 1:nrow(APSESMsaSet)) {
  seq <- as.character(APSESMsaSet@unmasked[i])
  msaMatrix[i, ] <- unlist(strsplit(seq, ""))
}

# inspect the result
msaMatrix[1:5, ]

# Now let's make a logical vector with an element
# for each column that selects which columns should
# be masked out.

# To count the number of elements in a vector, R has
# the table() function. For example ...
table(msaMatrix[ , 1])
table(msaMatrix[ , 10])
table(msaMatrix[ , 20])
table(msaMatrix[ , 30])


# Since the return value of table() is a named vector, where
# the name is the element that was counted in each slot,
# we can simply get the counts for hyphens from the
# return value of table(). We don't even need to assign
# the result to an intermediate variable, but we
# can attach the selection via square brackets,
# i.e.: ["-"],  directly to the function call:
table(msaMatrix[ , 1])["-"]

# ... to get the number of hyphens. And we can compare
# whether it is eg. > 4.
table(msaMatrix[ , 1])["-"] > 4

# Thus filling our logical vector is really simple:

# initialize the mask
colMask <- logical(lenAli)

# define the threshold for rejecting a column
limit <- round(nrow(APSESMsaSet) * (2/3))

# iterate over all columns, and write TRUE if there are less-or-equal to "limit"
# hyphens, FALSE if there are more.
for (i in 1:lenAli) {
  count <- table(msaMatrix[ , i])["-"]
  if (is.na(count)) { # No hyphen
    count <- 0
  }
  colMask[i] <- count <= limit
}

# inspect the mask
colMask

# How many positions were masked? R has a simple trick
# to count the number of TRUE and FALSE in a logical
# vector. If a logical TRUE or FALSE is converted into
# a number, it becomes 1 or 0 respectively. If we use
# the sum() function on the vector, the conversion is
# done implicitly. Thus ...
sum(colMask)

# ... gives the number of TRUE elements.

cat(sprintf("We are masking %4.2f %% of alignment columns.\n",
            100 * (1 - (sum(colMask) / length(colMask)))))


# Next, we use colMask to remove the masked columns from the matrix
# in one step:
maskedMatrix <- msaMatrix[ , colMask]

# check:
ncol(maskedMatrix)


# ... then collapse each row back into a sequence ...

apsMaskedSeq <- character()
for (i in 1:nrow(maskedMatrix)) {
  apsMaskedSeq[i] <- paste(maskedMatrix[i, ], collapse="")
}
names(apsMaskedSeq) <- rownames(maskedMatrix)

# ... and read it back into an AAStringSet object

apsMaskedSet <- AAStringSet(apsMaskedSeq)

# inspect ...
writeSeqSet(apsMaskedSet, format = "ali")



# Step 2. Turn this code into a function...

# Even though the procedure is simple, doing this more than once is tedious and
# prone to errors. I have assembled the steps we just went through into a
# function maskSet() and put it into the utilities.R file, from where it has
# been loaded when you started this sesssion.

maskSet

# Check that the function gives identical results
# to what we did before by hand:
identical(apsMaskedSet, maskSet(APSESMsaSet))

# The result must be TRUE. If it's not TRUE you have
# an error somewhere.

# We save the aligned, masked domains to a file in multi-FASTA format.
writeSeqSet(maskSet(APSESMsaSet), file = "APSES.mfa",   format = "mfa")



# = 1 Tasks




# [END]
