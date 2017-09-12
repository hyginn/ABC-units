# BIN-ALI-MSA.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-MSA unit.
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

# = 1 Multiple sequence alignment

# We will compute a multiple sequence alignment using the "muscle" algorithm
# which is available throught the Bioconductor msa package.

if (!require(Biostrings, quietly=TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
}
data(BLOSUM62)

if (!require(msa, quietly=TRUE)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("msa")
  library(msa)
}

# If the installation asks you if you want to update older packages, always
# answer "a" for "all" unless you have an important reason not to. But if the
# installer asks you whether you want to compile from source, answer"n" for "no"
# unless you need new functionality in a particular bleeding-edge version of a
# package.

help(package = "msa")

# We have used biostrings' AAString() function before; for multiple
# alignments we need an AAStringSet(). We can simply feed it
# a vector of sequences:

# Let's make a shorthand for selection of Mbp1 proteins from our database: a
# vector of indices for all of the rows in the protein table that contain
# proteins whose name begins with MBP1.
iMbp1Proteins <- grep("^MBP1_", myDB$protein$name)

# Next we create an AAStringSet for all of those proteins
seqSet <- AAStringSet(myDB$protein$sequence[iMbp1Proteins])

# ... and align (which is very simple) ...
msaMuscle(
  seqSet,
  order = "aligned")

# ... but to help us make sense of the alignment we need to add the names for
# the sequences. names for a seqSet object are held in the ranges slot...

seqSet@ranges@NAMES <- myDB$protein$name[iMbp1Proteins]

seqSet

# This little step of adding names is actually really
# very important. That's because the aligned sequences
# are meaningless strings of characters unless we can
# easily identify their biological relationships.
# Creating MSAs that are only identified by e.g. their
# RefSeq ids is a type of cargo-cult bioinformatics
# that we encounter a lot. The point of the alignment
# is not to create it, but to interpret it!


# Let's align again, and assign the result ...
msa1 <-  msaMuscle(
  seqSet,
  order = "aligned")

msa1

# ... or to see the whole thing (cf. ?MsaAAMultipleAlignment ... print method):
print(msa1, show=c("alignment", "complete"), showConsensus=FALSE)


# You see that the alignment object has sequence strings with hyphens as
# indel-characters. The names are printed to the console. And you also see that
# the order has not been preserved, but the most similar sequences are now
# adjacent to each other.

# Lets write the alignment to one of the common file formats: a multi-fasta
# file.

# Why oh why does the msa package not have a function to do this !!!
# Like, seriously ...

# ==== writeMFA

# Here is our own function to write a AAStringSet object to a multi-FASTA file.
writeMFA <- function(ali, file, blockSize = 50) {
  if (missing(ali)) {
    stop("Input object missing from arguments with no default.")
  }
  if (missing(file)) {
    writeToFile <- FALSE
  }
  else {
    writeToFile <- TRUE
    sink(file) # divert output to file
  }
  # Extract the raw data from the objects depending on
  # their respective class and put this
  # into a named vector of strings.
  if (class(ali)[1] == "MsaAAMultipleAlignment") {
    strings <- character(nrow(ali))
    for (i in 1:nrow(ali)) {
      strings[i] <- as.character(ali@unmasked[i])
      names(strings)[i] <- ali@unmasked@ranges@NAMES[i]
    }
  }
  else if (class(ali)[1] == "AAStringSet") {
    strings <- character(length(ali))
    for (i in 1:length(ali)) {
      strings[i] <- as.character(ali[i])
      names(strings)[i] <- ali@ranges@NAMES[i]
    }
  }
  else {
    stop(paste("Input object of class",
               class(ali)[1],
               "can't be handled by this function."))
  }


  for (i in 1:length(strings)) {
    # output FASTA header
    cat(paste(">",
              names(strings)[i],
              "\n",
              sep=""))
    # output the sequence block by block ...
    nLine <- 1
    from <- 1
    while (from < nchar(strings[i])) {
      to <- from + blockSize - 1
      cat(paste(substr(strings[i], from, to), "\n", sep=""))
      from <- to + 1
    }
    cat("\n") # output an empty line
  }
  if (writeToFile) {
    sink()  # Done. Close the diversion.
  }
}

# confirm that the function works
writeMFA(seqSet)
writeMFA(msa1)

# We could use this function to write the raw and aligned sequences to file like
# so:

# writeMFA(seqSet, file = "APSES_proteins.mfa")
# writeMFA(msa1, file = "APSES_proteins_muscle.mfa")

# ...  but since we don't actually need the data on file now, just copy the code
# that would create the msa to your myCode file so you can quickly reproduce it.

# == Task:
# Print the output of print(msa1) on a sheet of paper and bring it to
# class for the next quiz.

# That'll do.



# = 1 Tasks




# [END]
