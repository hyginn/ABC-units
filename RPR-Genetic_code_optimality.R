# tocID <- "RPR-Genetic_code_optimality.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Genetic_code_optimality unit.
#
# Version:  1.4
#
# Date:     2017-10  -  2022-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.4    2022: rewrite for a three-phase approach, more principled
#                        separation of concerns, and more principled
#                        quantification of similarity (with sections
#                        adapted from CSB195 code)
#           1.3    2020 Maintenance
#           1.2    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.1      Update set.seed() usage
#           1.0.1    Fixed two bugs discovered by Suan Chin Yeo.
#           1.0      New material.
#
#
# TODO:
#       ...
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
#TOC>   1        Context                                          57
#TOC>   1.1        Designing a computational experiment          100
#TOC>   1.2        Designing code                                114
#TOC>   2        First Phase: concept                            136
#TOC>   3        Setting up the tools                            151
#TOC>   3.1        Natural and alternative genetic codes         154
#TOC>   3.2        Effect of mutations                           236
#TOC>   3.2.1          Randomly mutate                           254
#TOC>   3.2.2          Translate                                 282
#TOC>   3.2.3          quantify mutation effects                 300
#TOC>   3.3        Computing pairwise similarity                 343
#TOC>   4        Run the experiment                              374
#TOC> 
#TOC> ==========================================================================


# =    1  Context  =============================================================
#
#
#   +-------------------------------------------+
#   |                                           |
#   |  Is the genetic code the best it can be?  |
#   |                                           |
#   +-------------------------------------------+
#
# This is actually an ill-posed question, as long as we do not have a precise
# definition of what we mean by "the best". But we can actually come up with
# a definition in the following way:
#
#   1. Obviously, the genetic code is non-random
#   2. If it is not random, it has evolved under selective pressure.
#   3. We can construct computational experiments about selective pressure. We
#      define a particular type of selective pressure, evaluate how well
#      the natural code does under that criterion, and compare its performance
#      to synthetic, simulated codes. If it turns out the natural code is much
#      more efficient than our simulated alternatives, that selective pressure
#      becomes a candidate hypothesis on how the code has actually evolved.
#   4. The question then becomes:
#        Does the genetic code respond to [the selective pressure we have
#        just defined] better than [almost] all alternatives?
#      ... and this question can be decided by a computational experiment!

# Here are some types of selective pressure we might formulate. If you can
# come up with additional ideas, share them, so we can add them to the list.
#  -  Conservative: minimize the effect of mutations.
#  -  Explorative: maximize the chance of a mutation being useful
#  -  Evolvable: the code must be derivable from simpler precursors
#  -  Efficient: the code should store the greatest amount of information in the
#  -             shortest DNA sequence
#
# Obviously, these criteria are not entirely orthogonal. In this code unit
# we will only explore the first one: the universal genetic code is conservative
# and minimizes the effect of mutations in a sequence.

# What we do here is deliberately simple - you are encouraged to think of
# refinements; this code is mainly intended to demonstrate some R programming
# techniques, in preparation for later units.


# ==   1.1  Designing a computational experiment  ==============================

# Computational experiments are conducted like wet-lab experiments. We begin
# with a hypothesis, then define the observables that relate to the hypothesis,
# then define the measures we apply to observations, and finally we interpret
# our observations. If we want to test the hypothesis that the genetic code
# has evolved to be conservative ...

#  - we define the observables: the effect of mutations in a given
#      sequence, given the natural and possible alternative codes;
#  - we define a quantitative measure of the effect of mutations;
#  - then we compute alternative codes and interpret the results.


# ==   1.2  Designing code  ====================================================
#
# This course has a separate section on software design:
#   http://steipe.biochemistry.utoronto.ca/bio/FND-CSC-Software_development.html
# ... but in principle we always approach the process in three phases:
#
# - First phase:  Write down the skeleton of the process in pseudocode
# - Second phase: For each step of the pseudocode, define the data structures
#                 and functions it requires. The functions might just be mock up
#                 functions initially that merely accept input and return
#                 dummy output.
# - Third phase:  For each step of the pseudocode, implement the actual process.
#                 As you do this, the pseudocode may need to be refined, the
#                 functions are implemented, detail added to the datastructures
#                 etc.
# Iterate.
#
# Naturally, the *first phase* section of your document will turn into
# documentation of concept and intent, the *second phase* section will define
# your functions, and the *third phase* will implement your actual process.


# =    2  First Phase: concept  ================================================
#
# We follow the ideas above:

#   Represent the universal genetic code (UGC) as a data structure
#   for N trials
#     Construct a random genetic code (RGC)
#     for all single point mutations in a given sequence
#       - compute a similarity value for the two amino acids that are encoded
#       - add it to a total (smaller sums are more similar codes)
#     store each sum
#   Tabulate distribution of sums for all simulated RGC
#   Report where the UGC fall on this distribution


# =    3  Setting up the tools  ================================================


# ==   3.1  Natural and alternative genetic codes  =============================


#   Represent the universal genetic code (UGC) as a data structure
#   --------------------------------------------------------------
#
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

uGC <- Biostrings::GENETIC_CODE


#     Construct a random genetic code (RGC)
#     -------------------------------------

# There are many ways to generate random, alternative codes. "Fully" random
# would not guarantee that all amino acids are represented and that would be
# non-biological. Thus we would like to conserve the alphabet. We could also
# conserve redundancy of the code. Preserving redundancy of codons, may reflect
# some form of symmetry breaking that ignores the third nucleotide of a codon
# for the most part. We'll implement that as well.

randomGC <- function(GC) {
  # Return a genetic code with randomly assigned amino acids.
  # Parameters:
  #    GC   named chr  length-64 character vector of 20 amino acid one-letter
  #                       codes plus "*" (stop), named with the codon triplet.
  # Value:  named chr  same vector with random amino acid assignments in which
  #                       every amino acid and "*" is encoded at least once.

  aa <- unique(GC)                           # the amino acids in the input code
  GC[1:64] <- sample(aa, 64, replace=TRUE)   # random code
  while(length(unique(GC)) < length(aa)) {   # We could end up with a code that
                                             # does not contain all amino acids,
                                             # then we sample() again.
    GC[1:64] <- sample(aa, 64, replace = TRUE)
  }
  return(GC)
}

isoredGC <- function(GC) {
  # Return a genetic code with the same redundancy as the original code -
  # "iso-redundant".
  # Parameters:
  #    GC   named chr  length-64 character vector of 20 amino acid one-letter
  #                       codes plus "*" (stop), named with the codon triplet.
  # Value:  named chr  same vector with random amino acid assignments where the
  #                       amino acids have been swapped and permuted.

  codonNames <- names(GC)

  # 1: use chartr() to replace each amino acid with a different letter
  alphabet <- unique(GC)
  alphrand <- sample(alphabet)
  irGC <- chartr(paste(alphabet, collapse = ""),
                 paste(alphrand, collapse = ""),
                 GC)
  # 2: permute the ordering of the amino acids
  irGC <- sample(irGC)
  # 3: restore the original order of codons
  names(irGC) <- codonNames

  return(irGC)
}

# for consistency, we'll also define an identity function - it merely
# return the original - i.e. it does nothing.

idGC <- function(GC) {
  return(GC)
}


# ==   3.2  Effect of mutations  ===============================================


# To evaluate the effects of mutations we will do the following:
#   - we take a nucleotide sequence (Mbp1 will do just nicely);
#   - we mutate it randomly;
#   - we translate it back to amino acids;
#   - we count the number of mutations and evaluate their severity.

# Mbp1 nucleotide sequence taken from https://www.ncbi.nlm.nih.gov/gene/851503
# (click on FASTA to get the correct range of 352877-355378 on chromosome IV,
#  copy, paste into a text object, sanitize, strsplit(), and saveRDS() ...)

myDNA <- readRDS("./data/MBP1SACCE_gene.rds")
head(myDNA)
tail(myDNA)
length(myDNA) / 3   # 833 amino acids + stop codon

# ===   3.2.1  Randomly mutate                      

# To mutate, we randomly choose a position, and randomly replace it with one
# of the three other nucleotides.

randMut <- function(myC) {
  # Parameter:
  # myC     chr     a codon - a vector of three characters from {A C G T}
  #                 no checking for correct case or alphabet is done.
  # Value:  chr     a codon with a single point mutation relative to myC

  nuc <- c("A", "C", "G", "T")     # define the alphabet

  mutC <- myC

  iMut <- sample(1:length(myC), 1)           # choose one position to mutate
  mutNuc <- sample(nuc[nuc != myC[iMut]], 1) # chose a different nucleotide
  mutC[iMut] <- mutNuc                        # replace the original

  return(mutC)
}

# Verify  (Hm. How does this even work?)
cat(sprintf("codon: %s mutated: %s\n",
    paste(x <- sample(c("A","C","G","T"), 3, replace=TRUE), collapse=""),
    paste(randMut(x), collapse = "")))


# ===   3.2.2  Translate                            

tra <- function(myC, GC) {
  # translate a codon to an amino acid as defined in the genetic code GC
  # Parameters:
  #   myC  chr   a codon - a vector of three characters from {A C G T}
  #   GC   chr   a genetic code with all 64 permuations of {A C G T}
  #
  # Value:
  #        chr   an amino acid or stop as defined in GC

  cod <- paste(myC, collapse = "")
  myAA <- GC[cod]  # translate via name

  return(myAA)
}


# ===   3.2.3  quantify mutation effects            

# How do we evaluate the effect of the mutation? We'll take a simple approach:
# we take three vectors of physicochemical properties from the seqinr:: aaindex
# tables. This defines a point in a 3D feature space for every amino acid. Then
# we define "similarity" as the distance between two points.

# fetch three feature vectors from seqinr::aaindex
library(seqinr)
data(aaindex)

# 544: Hydrophobicity index (Fasman, 1989)
myF <- data.frame(f1=aaindex[[544]]$I)

# 515: Average volumes of residues (Pontius et al., 1996)
myF$f2 <- aaindex[[515]]$I

# 8: Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)
myF$f3 <- aaindex[[8]]$I  # flexibility

# To avoid problems with different absolute values of our features, we
# scale() to center the columns (each column's mean becomes 0) and scales
# the values so that each vector's standard deviation is one. This type of
# "normalization" is called "standardization". Scale can work on all columns
# of a matrix-like object at once.

myF <- as.data.frame(scale(myF))

# Change the rownames of myF to single-letter code
rownames(myF) <- seqinr::a(rownames(myF))  # reassign them

# Finally: what should the similarity of an amino acid and the stop codon be?
# We can define this in many ways. For now, I just set the value to
# zero.
myF <- rbind(myF, c(0, 0, 0))

# add the stop-codon rowname
rownames(myF)[nrow(myF)] <- "*"

# Verify...
tail(myF)


# ==   3.3  Computing pairwise similarity  =====================================


sim <- function(a, b, F = myF) {
  # a, b:  amino acids in one-letter code
  # F:     an amino acid feature table (default: myF)
  # value: the Euclidian distance between the two vectors F[a] and F[b]

  d <- sqrt(sum( (F[a, ] - F[b, ])^2 ))  # distance in Euclidian space
  return(d)
}

# Verify:
sim("H", "H")   # identical
sim("K", "E")   # rather similar
sim("S", "W")   # quite dissimilar

# Of course we could take additional features into account, such as charge,
# or H-bond propensity - and perhaps add heuristics, such as: proline is always
# bad in secondary structure, charged amino acids are terrible in the folded
# core of a protein, replacing a small by a large amino acid in the core is very
# disruptive ... etc.

# However, we cannot use a mutation data matrix for this experiment! While
# empirical mutation probabilities are superbly suited to estimate evolutionary
# relationships, here we are trying to evaluate effects of random mutations on
# genetic codes, and that would make our reasoning circular - we would
# discover that the natural genetic code is optimal ... because it is most
# similar to the natural genetic code. That would be Cargo Cult bioinformatics.


# =    4  Run the experiment  ==================================================

# To run the experiment, we build a small function that ties
# together the workflow we have decided on:

runMutations <- function(N,                # Number of trials
                         DNA,              # Sequence to mutate
                         GC,               # Genetic code to use
                         fSynth = idGC) {  # function to change the code

  # Value: a vector of length N containing the summed similarity values
  # for a fully mutated sequence


  idxCodons <- seq(1, length(DNA), by=3) # index vector of codon start positions
  results <- numeric(N)                  # vector to store results

  cat(sprintf("\nRunning %d trials...\n", N))

  for (i in 1:N) {
    pBar(i, N)                                  # show a progress bar

    fx <- numeric(length(idxCodons))            # initialize a vector to store
                                                # the effects of one trial

    thisGC <- fSynth(GC)                        # make a synthetic code

    for (j in 1:length(idxCodons)) {            # for each codon
      myPos <- idxCodons[j]:(idxCodons[j] + 2)  # determine position
      og <- DNA[myPos]                          # get original codon
      mu <- randMut(og)                         # mutate the codon
      fx[j] <- sim(tra(og, thisGC),             # compute and store the effect
                   tra(mu, thisGC))
    }
    results[i] <- sum(fx)                       # save the sum of all effects
                                                # for each trial
  }
  return(results)
}

# Make sure the values are correctly defined:

# Fetch the standard Genetic code from Biostrings::
stdCode <- Biostrings::GENETIC_CODE

# Fetch the nucleotide sequence for MBP1:
myDNA <- readRDS("./data/MBP1SACCE_gene.rds")

# The amino acid feature table myF has been defined above. We won't redefine it.
#

N <- 200
fxStdCode <- runMutations(N, DNA = myDNA, GC = stdCode)

# Note: this is really quite slow. The reason is we are doing every single
# mutation on its own - all N * 834 of them - and not taking advantage of
# vectorized functions and precomputed datastructures at all. But coding it more
# efficiently would obscure the principle.

myBreaks <- seq(1000, 2100, by = 30)
hist(fxStdCode,
     breaks = myBreaks,
     col = "#8bbccc",
     xlim = c(1000, 2200),
     ylim = c(0, 100),
     main = "Standard vs. Synthetic Genetic Codes",
     xlab = "Summed similarities",
     ylab = sprintf("N observations (of %s)", N))


# Compare this to random genetic codes: fully random
#fxRndCode <- runMutations(N, DNA = myDNA, GC = stdCode, fSynth = randomGC)

hist(fxRndCode,
     breaks = myBreaks,
     col = "#4c679388",
     add = TRUE)

# ... and iso-redundant
#fxIsoCode <- runMutations(N, DNA = myDNA, GC = stdCode, fSynth = isoredGC)

hist(fxIsoCode,
     col = "#e0afe388",
     breaks = myBreaks,
     add = TRUE)


# Interpret.

# These are simple experiments, under assumptions that can be refined in
# meaningful ways. Yet, even those simple computational experiments show
# that the Standard Genetic Code has features that one would predict if
# it has evolved under selective pressure to minimize the effects of mutations.
# Gradual change under mutation is beneficial to evolution, disruptive
# change is not.


# [END]
