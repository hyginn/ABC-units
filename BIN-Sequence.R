# BIN-Sequence.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-Sequence unit.
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
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...
#
# ==============================================================================

#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                          Line
#TOC> ----------------------------------------------
#TOC>   1        Prepare                          47
#TOC>   2        Storing Sequence                 61
#TOC>   3        String properties                90
#TOC>   4        Substrings                       97
#TOC>   5        Creating strings: sprintf()     103
#TOC>   6        Changing strings                134
#TOC>   6.1      stringi and stringr             162
#TOC>   6.2      dbSanitizeSequence()            172
#TOC>   7        Permuting and sampling          184
#TOC>   7.1      Permutations                    191
#TOC>   7.2      Sampling                        234
#TOC>   7.2.1    Equiprobable characters         236
#TOC>   7.2.2    Defined probability vector      271
#TOC>   8        Tasks                           299
#TOC>
#TOC> ==========================================================================


# =    1  Prepare  =============================================================

# Much basic sequence handling is supported by the Bioconductor package
# Biostrings.

if (! require(Biostrings)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biostrings")
  library(Biostrings)
}


# =    2  Storing Sequence  ====================================================


# Sequences can be represented and stored as vectors of single characters ...
(v <- c("D", "I", "V", "M", "T", "Q"))

# ... as strings ...
(s <- "DIVMTQ")

# ... or as more complex objects with rich metadata e.g. as a Biostrings
# DNAstring, RNAstring, AAString, etc.
(a <- AAString("DIVMTQ"))

# ... and all of these representations can be interconverted:

# string to vector ...
unlist(strsplit(s, ""))

# vector to string ...
paste(v, sep = "", collapse = "")

# ... and AAstring to plain string.
as.character(a)

# Since operations with character vectors trivially follow all other vector
# conventions and syntax, and we will look at Biostrings methods in more
# detail in a later unit, we will focus on basic strings in the following.


# =    3  String properties  ===================================================


length(s) # why ???
nchar(s)  # aha


# =    4  Substrings  ==========================================================

# Use the substr() function
substr(s, 2, 4)


# =    5  Creating strings: sprintf()  =========================================


# Sprintf is a _very smart, very powerful function and has cognates in all
# other programming languages. It has a small learning curve, but it's
# totally worth it:
# the function takes a format string, and a list of other arguments. It returns
# a formatted string. Here are some examples - watch carefully for sprintf()
# calls in other code.

sprintf("Just a string.")
sprintf("A string and the number %d.", 5)
sprintf("More numbers: %d ate %d.", 7, 9) # Sorry
sprintf("Pi is ~ %1.2f ...", pi)
sprintf("or more accurately ~ %1.11f.", pi)
x <- "bottles of beer"
n <- 99
sprintf("%d %s on the wall, %d %s - \ntake %s: %d %s on the wall.",
        n, x, n, x, "one down, and pass it around", n-1, x)

# Note that in the last example, the value of the string was displayed with
# R's usual print-formatting function and therefore the line-break "\n" did
# not actually break the line. To have line breaks, tabs etc, you need to use
# cat() to display the string:

for (i in 99:95) {
  cat(sprintf("%d %s on the wall, %d %s - \ntake %s: %d %s on the wall.\n\n",
              i, x, i, x, "one down, and pass it around", i-1, x))
}


# =    6  Changing strings  ====================================================

# Changing case
tolower(s)
toupper(tolower(s))

#reverse
reverse(s)

# substituing characters
(s <- gsub("IV", "i-v", s))  # gsub can change length, first argument is
                             # a "regular expression"!

# I use it often to delete characters I don't want ...
(s <- gsub("-", "", s))

# For example: clean up a sequence
# copy/paste from UniProt
(s <- "        10         20         30         40         50
MSNQIYSARY SGVDVYEFIH STGSIMKRKK DDWVNATHIL KAANFAKAKR ")


# remove numbers
(s <- gsub("[0-9]", "", s))

# remove "whitespace" (spaces, tabs, line breaks)...
(s <- gsub("\\s", "", s))

# ==   6.1  stringi and stringr  ===============================================

# But there are also specialized functions eg. to remove leading/trailing
# whitespace which may be important to sanitize user input etc. Have a look at
# the function descriptions for the stringr and the stringi package. stringr is
# part of the tidyverse, and for the most part a wrapper for stringi functions.
# https://github.com/tidyverse/stringr



# ==   6.2  dbSanitizeSequence()  ==============================================

# In our learning units, we use a function dbSanitizeSequence() to clean up
# sequences that may be copy/pasted from Web-sources

s <- ">FASTA header will be removed
10         20         30         40         50
MSNQIYSARY SGVDVYEFIH STGSIMKRKK DDWVNATHIL KAANFAKAKR "

dbSanitizeSequence(s)


# =    7  Permuting and sampling  ==============================================


# An important aspect of working with strings is generating random strings
# with given statistical properties: reference items to evaluate significance.


# ==   7.1  Permutations  ======================================================


# One way to produce such reference items is to permute a string. A permuted
# string has the same composition as the original, but all positional
# information is lost. The sample() function can be used to permute:

# This is the sequence of the ompA secretion signal
(s <- unlist(strsplit("MKKTAIAVALAGFATVAQA", "")))

(x <- sample(s, length(s)))  # permuted

# Here's a small example how such permuted strings may be useful. As you look
# at the ompA sequence, you suspect that the two lysines near the +-charged
# N-terminus may not be accidental, but selected for a positively charged
# N-terminus. What is the chance that such a sequence has two lysines close to
# the N-terminus simply by chance? Or put differently: what is the average
# distance of two lysines in such a sequence to the N-terminus. First, we
# need an expression that measures the distance. A simple use of the which()
# function will do just fine.

which(s == "K")        # shows they are in position 2 and 3, so ...
mean(which(s == "K"))  # ... gives us the average, and ...
mean(which(x == "K"))  # ... gives us the average of the permuted sequence.

# So what does the distribution look like? Lets do 10,000 trials.

(s <- unlist(strsplit("MKKTAIAVALAGFATVAQA", "")))
N <- 10000
d <- numeric(N)
set.seed(112358)
for (i in 1:N) {
  d[i] <- mean(which(sample(s, length(s)) == "K"))
}
hist(d, breaks = 20)
abline(v = 2.5, lwd = 2, col = "firebrick")
sum(d <= 2.5) # 276. 276 of our 10000 samples are just as bunched near the
              # N-terminus or more. That's just below the signifcance
              # threshold of 5 %. It's a trend, but to be sure we are looking
              # at a biological effect we would need to see more
              # sequences.


# ==   7.2  Sampling  ==========================================================

# ===  7.2.1  Equiprobable characters

# Assume you need a large random-nucleotide string for some statistical model.
# How to create such a string? sample() can easily create it:

nuc <- c("A", "C", "G", "T")
N <- 100
set.seed(16818)
v <- sample(nuc, N, replace = TRUE)
(mySeq <- paste(v, collapse = ""))

# What's the GC content?
table(v)
sum(table(v)[c("G", "C")]) # 51 is close to expected

# What's the number of CpG motifs? Easy to check with the stringi
# stri_match_all() function

if (! require(stringi)) {
  install.packages("stringi")
  library(stringi)
}

(x <- stri_match_all(mySeq, regex = "CG"))
length(unlist(x))

# Now you could compare that number with yeast DNA sequences, and determine
# whether there are more or less CpG motifs than expected by chance.
# (cf. https://en.wikipedia.org/wiki/CpG_site)
# But hold on: is that a fair comparison? sample() gives us all four nucleotides
# with the same probability. But the yeast genomic DNA GC content is only
# 38%. So you would expect fewer CpG motifs based on the statistical properties
# of the smaller number of Cs and Gs - before biology even comes into play. How
# do we account for that?

# ===  7.2.2  Defined probability vector

# This is where we need to know how to create samples with specific probability
# distributions. A crude hack would be to create a sampling source vector with
# 19 C, 19 G, 31 A and 31 T
c(rep("C", 19), rep("G", 19), rep(c("A"), 31), rep(c("T"), 31))
# ... but that doesn't scale if the numeric accuracy needs to be higher.
#
# However sample() has an argument that takes care of that: you can explicitly
# specify the probabilities with which each element of the the sampling vector
# should be chosen:

nuc <- c("A", "C", "G", "T")
N <- 100
set.seed(16818)
myProb <- c(0.31, 0.19, 0.19, 0.31)  # sampling probabilities
v <- sample(nuc, N, prob = myProb, replace = TRUE)
(mySeq <- paste(v, collapse = ""))

# What's the GC content?
table(v)
sum(table(v)[c("G", "C")]) # Close to expected

# What's the number of CpG motifs?
(x <- stri_match_all(mySeq, regex = "CG"))
# ... not a single one in this case.


# =    8  Tasks  ===============================================================

# Task: Phone numbers that are entered into Web forms can come in many
#         different formats. Write a function sanitizePhone() that accepts
#         a single object as input and returns a single string of only numbers.

          if (! require(testthat)) {
            install.packages("testthat")
            library(testthat)
          }

          sanitizePhone <- function(s) {
            # ... your function code here
          }

          # All tests must pass!
          s <- "1-858 651-5050"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "1 858 651 5050"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "+1 (858) 651-5050"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "18586515050"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "1 858 6515050"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "1.858.651.5050"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "1\t8 5 8\t6 5 1-5 0 5 0"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "1n8e5v8e6r5 1g5o0n5n0a"
          expect_equal(sanitizePhone(s), "18586515050")
          s <- "IDK"
          expect_equal(sanitizePhone(s), "")
          s <- ""
          expect_equal(sanitizePhone(s), "")
          s <- pi
          expect_equal(sanitizePhone(s), "314159265358979")



# [END]