# tocID <- "RPR-Genetic_code_optimality.R"
#
# ---------------------------------------------------------------------------- #
#  PATIENCE  ...                                                               #
#    Do not yet work wih this code. Updates in progress. Thank you.            #
#    boris.steipe@utoronto.ca                                                  #
# ---------------------------------------------------------------------------- #
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Genetic_code_optimality unit.
#
# Version:  1.2
#
# Date:     2017  10  -  2019  01
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.1      Update set.seed() usage
#           1.0.1    Fixed two bugs discovered by Suan Chin Yeo.
#           1.0      New material.
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
#TOC>   Section  Title                                          Line
#TOC> --------------------------------------------------------------
#TOC>   1        Designing a computational experiment             57
#TOC>   2        Setting up the tools                             73
#TOC>   2.1        Natural and alternative genetic codes          76
#TOC>   2.2        Effect of mutations                           134
#TOC>   2.2.1          reverse-translate                         145
#TOC>   2.2.2          Randomly mutate                           170
#TOC>   2.2.3          Forward- translate                        195
#TOC>   2.2.4          measure effect                            213
#TOC>   3        Run the experiment                              260
#TOC>   4        Task solutions                                  356
#TOC>
#TOC> ==========================================================================


# This unit demonstrates R code to simulate alternate genetic codes and evaluate
# their robsustness to code changes. The approaches are quite simple and you
# will be able to come up with obvious refinements; the point of this code is to
# demonstrate some R programming techniques, in preparation for more
# sophisticated questions later.


# =    1  Designing a computational experiment  ================================

# Computational experiments are conducted like wet-lab experiments. We begin
# with a hypothesis, then define the observables that relate to the hypothesis,
# then define the measures we apply to observations, and finally we interpret
# our observations. If we want to learn something about the evolution of the
# genetic code ...

#  - we construct a hypothesis such as: the genetic code has evolved so as to
#      minimize the effect of mutations;
#  - we define the observables: the effect of mutations in
#      sequences, given the natural and possible alternative codes;
#  - we define the measures to quantify the effect of mutations;
#  - then we compute alternatives and interpret the results.


# =    2  Setting up the tools  ================================================


# ==   2.1  Natural and alternative genetic codes  =============================

# Load genetic code tables from the Biostrings package
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


# There are many ways to generate alternative codes. The simplest way is to
# randomly assign amino acids to codons. A more sophisticated way is to keep the
# redundancy of codons intact, since it may reflect some form of symmetry
# breaking that ignores the third nucleotide of a codon for the most part;
# therefore we only replace the amino acids of the existing code with random
# others. Here are two functions that implement these two ideas about alternate
# codes.

randomGC <- function(GC) {
  # Return a genetic code with randomly assigned amino acids.
  # Parameters:
  #    GC   named chr  length-64 character vector of 20 amino acid one-letter
  #                       codes plus "*" (stop), named with the codon triplet.
  # Value:  named chr  same vector with random amino acid assignments in which
  #                       every amino acid and "*" is encoded at least once.

  aa <- unique(GC)                           # the amino acids in the input code
  GC[1:64] <- sample(aa, 64, replace = TRUE) # random code
  while(length(unique(GC)) < length(aa)) {   # We could end up with a code that
                                             # does not contain all amino acids,
                                             # then we sample() again.
    GC[1:64] <- sample(aa, 64, replace = TRUE)
  }
  return(GC)
}

swappedGC <- function(GC) {
  # Return a genetic code with randomly swapped amino acids.
  # Parameters:
  #    GC   named chr  length-64 character vector of 20 amino acid one-letter
  #                       codes plus "*" (stop), named with the codon triplet.
  # Value:  named chr  same vector with random amino acid assignments where the
  #                       amino acids have been swapped.

  aaOrig <- unique(GC)                       # the amino acids in the input code
  aaSwap <- sample(aaOrig, length(aaOrig))   # shuffled
  names(aaSwap) <- aaOrig                    # name them after the original
  GC[1:64] <- aaSwap[GC]                     # replace original with shuffled

  return(GC)
}


# ==   2.2  Effect of mutations  ===============================================


# To evaluate the effects of mutations we will do the following:
#   - we take an amino acid sequence (Mbp1 will do just nicely);
#   - we reverse-translate it into a nucleotide sequence;
#   - we mutate it randomly;
#   - we translate it back to amino acids;
#   - we count the number of mutations and evaluate their severity.


# ===   2.2.1  reverse-translate

# To reverse-translate an amino acid vector, we randomly pick one of its
# codons from a genetic code, and assemble all codons to a sequence.

traRev <- function(s, GC) {
  # Parameters:
  #      s   chr   a sequence vector
  #      GC  chr   a genetic code
  # Value:
  #      A reverse-translated vector of codons
  vC <- character(length(s))

  for (i in seq_along(s)) {
    codon <- names(GC)[GC == s[i]]   # get all codons for this AA
    if (length(codon) > 1) {         # if there's more than one ...
      codon <- sample(codon, 1)      # pick one at random ...
    }
    vC[i] <- codon                   # store it
  }

  return(vC)
}


# ===   2.2.2  Randomly mutate

# To mutate, we split a codon into it's three nucleotides, then randomly replace
# one of the three with another nucleotide.

randMut <- function(vC) {
  # Parameter:
  #    vC   chr     a vector of codons
  # Value:  chr     a vector of codons with a single point mutation from vC

  nuc <- c("A", "C", "G", "T")

  for (i in seq_along(vC)) {
    triplet <- unlist(strsplit(vC[i], ""))         # split into three nucl.
    iNuc <- sample(1:3, 1)                         # choose one of the three
    mutNuc <- sample(nuc[nuc != triplet[iNuc]], 1) # chose a mutated nucleotide
    triplet[iNuc] <- mutNuc                        # replace the original
    vC[i] <- paste0(triplet, collapse = "")        # collapse it to a codon
  }
  return(vC)

}



# ===   2.2.3  Forward- translate

traFor <- function(vC, GC) {
  # Parameters:
  #      vC   chr   a codon vector
  #      GC   chr   a genetic code
  # Value:
  #      A vector of amino acids
  vAA <- character(length(vC))

  for (i in seq_along(vC)) {
    vAA[i] <- GC[vC[i]]         # translate and store
  }
  return(vAA)

  }


# ===   2.2.4  measure effect

# How do we evaluate the effect of the mutation? We'll take a simple ad hoc
# approach: we divide amino acids into hydrophobic, hydrophilic, and neutral
# categories, according to their free energy of transfer from water to octanol:
aaHphobic <- c("M", "I", "L", "C", "W", "Y", "F")
aaHphilic <- c("E", "D", "Q", "N", "P", "K", "R")
aaNeutral <- c("A", "H", "T", "S", "V", "G")

# Then we will penalize as follows:
# Changes within one category: 0.1
# Changes from hydrophobic or hydrophilic to neutral or back: 0.3
# Changes from hydrophobic to hydrophilic or back: 1.0
# Changes to stop-codon: 3.0

evalMut <- function(nat, mut) {
  # Evaluate severity of mutations between amino acid sequence vectors nat and
  # mut in an ad hoc approach based on hydrophobicity changes.
  aaHphobic <- c("M", "I", "L", "C", "W", "Y", "F")
  aaHphilic <- c("E", "D", "Q", "N", "P", "K", "R")
  aaNeutral <- c("A", "H", "T", "S", "V", "G")

  penalties <- numeric(length(nat))
  lMut <- nat != mut    # logical TRUE for all mutated positions

  penalties[lMut & (nat %in% aaHphobic) & (mut %in% aaHphobic)] <- 0.1
  penalties[lMut & (nat %in% aaHphobic) & (mut %in% aaHphilic)] <- 1.0
  penalties[lMut & (nat %in% aaHphobic) & (mut %in% aaNeutral)] <- 0.3

  penalties[lMut & (nat %in% aaHphilic) & (mut %in% aaHphobic)] <- 1.0
  penalties[lMut & (nat %in% aaHphilic) & (mut %in% aaHphilic)] <- 0.1
  penalties[lMut & (nat %in% aaHphilic) & (mut %in% aaNeutral)] <- 0.3

  penalties[lMut & (nat %in% aaNeutral) & (mut %in% aaHphobic)] <- 0.3
  penalties[lMut & (nat %in% aaNeutral) & (mut %in% aaHphilic)] <- 0.3
  penalties[lMut & (nat %in% aaNeutral) & (mut %in% aaNeutral)] <- 0.1

  return(sum(penalties))
}

# A more sophisticated approach could take additional quantities into account,
# such as charge, size, or flexibility - and it could add heuristics, such as:
# proline is always bad in secondary structure, charged amino acids are terrible
# in the folded core of a protein, replacing a small by a large amino acid in
# the core is very disruptive ... etc.


# =    3  Run the experiment  ==================================================

# Fetch the standard Genetic code from Biostrings::

stdCode <- Biostrings::GENETIC_CODE

# Fetch the nucleotide sequence for MBP1:

myDNA <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")[-1]
myDNA <- paste0(myDNA, collapse = "")
myDNA <- as.character(Biostrings::codons(Biostrings::DNAString(myDNA)))
myDNA <- myDNA[-length(myDNA)]  # drop the stop codon

myAA <- traFor(myDNA, stdCode)

# Mutate and evaluate
set.seed(112358)
x <- randMut(myDNA)
set.seed(NULL)
x <- traFor(x, stdCode)
evalMut(myAA, x)  # 166.4

# Try this 200 times, and see how the values are distributed.
N <- 200
valSTDC <- numeric(N)

set.seed(112358)                   # set RNG seed for repeatable randomness
for (i in 1:N) {
  x <- randMut(myDNA)              # mutate
  x <- traFor(x, stdCode)     # translate
  valSTDC[i] <- evalMut(myAA, x)    # evaluate
}
set.seed(NULL)                     # reset the RNG

hist(valSTDC,
     breaks = 15,
     col = "palegoldenrod",
     xlim = c(0, 400),
     ylim = c(0, N/4),
     main = "Standard vs. Synthetic Genetic Code",
     xlab = "Mutation penalty")

# This looks like a normal distribution. Let's assume the effect of mutations
# under the standard genetic code is the mean of this distribution:
effectSTDC <- mean(valSTDC)  # 178.1

# Now we can look at the effects of alternate genetic codes:

set.seed(112358)
# choose a new code
GC <- randomGC(stdCode)
set.seed(NULL)

# reverse translate hypothetical sequence according to the new code
x <- traRev(myAA, GC)

x <- randMut(x)        # randomly mutate hypothetical nucleotide sequence
x <- traFor(x, GC)     # translate back, with the new code
evalMut(myAA, x)       # evaluate mutation effects: 298.5

# That seems a fair bit higher than what we saw as "effectUGC"
# Let's try with different genetic codes. 200 trials - but this time every trial
# is with a different, synthetic genetic code.

N <- 200
valXGC <- numeric(N)

set.seed(1414214)                # set RNG seed for repeatable randomness
for (i in 1:N) {
  GC <- randomGC(stdCode)   # Choose code
  x <- traRev(myAA, GC)          # reverse translate
  x <- randMut(x)                # mutate
  x <- traFor(x, GC)             # translate
  valXGC[i] <- evalMut(myAA, x)  # evaluate
}
set.seed(NULL)                   # reset the RNG

hist(valXGC,
     col = "plum",
     breaks = 15,
     add = TRUE)

# These two distributions are very widely separated!

# Task: Perform the same experiment with the swapped genetic code.
#       Compare the distributions. Interpret the result.


# These are simple experiments, under assumptions that can be refined in
# meaningful ways. Yet, even those simple computational experiments show
# that the Universal Genetic Code has features that one would predict if
# it has evolved under selective pressure to minimize the effects of mutations.
# Gradual change under mutation is benificial to evolution, disruptive
# change is not.


# =    4  Task solutions  ======================================================

N <- 200
valSGC <- numeric(N)

set.seed(2718282)                # set RNG seed for repeatable randomness
for (i in 1:N) {
  GC <- swappedGC(stdCode)  # Choose code
  x <- traRev(myAA, GC)          # reverse translate
  x <- randMut(x)                # mutate
  x <- traFor(x, GC)             # translate
  valSGC[i] <- evalMut(myAA, x)  # evaluate
}
set.seed(NULL)                   # reset the RNG

hist(valSGC,
     col = "#6688FF88",
     breaks = 15,
     add = TRUE)



# [END]
