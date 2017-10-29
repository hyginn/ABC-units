# BIN-ALI-Similarity.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-Similarity unit.
#
# Version:  1.0
#
# Date:     2017  10  20
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    Refactored for 2017; add aaindex, ternary plot.
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
#TOC>   Section  Title                    Line
#TOC> ----------------------------------------
#TOC>   1        Amino Acid Properties      43
#TOC>   2        Mutation Data matrix      163
#TOC>   3        Background score          205
#TOC> 
#TOC> ==========================================================================


# =    1  Amino Acid Properties  ===============================================

# A large collection of amino acid property tables is available via the seqinr
# package:

if (!require(seqinr)) {
  install.packages("seqinr")
  library(seqinr)
}
# Package information:
#  library(help = seqinr)       # basic information
#  browseVignettes("seqinr")    # available vignettes
#  data(package = "seqinr")     # available datasets

# A true Labor of Love has gone into the compilation of the seqinr "aaindex"
#  data:

?aaindex
data(aaindex)  # load the aaindex list from the package

length(aaindex)

# Here are all the index descriptions
for (i in 1:length(aaindex)) {
  cat(paste(i, ": ", aaindex[[i]]$D, "\n", sep=""))
}

# It's a bit cumbersome to search through the descriptions ... here is a
# function to make this easier:

searchAAindex <- function(patt) {
  # Searches the aaindex descriptions for regular expression "patt"
  # and prints index number and description.
  hits <- which(sapply(aaindex, function(x) length(grep(patt, x$D)) > 0))
  for (i in seq_along(hits)) {
    cat(sprintf("%3d\t%s\n", hits[i], aaindex[[ hits[i] ]]$D))
  }
}


searchAAindex("free energy")          # Search for "free energy"
searchAAindex("(size)|(volume)")      # Search for "size" or "volume":




# Let's examine ...
# ... a hydrophobicity index
(Y <- aaindex[[528]][c("D", "I")])

# ... a volume index
(V <- aaindex[[150]][c("D", "I")])

# ... and one of our own: side-chain pK values as reported by
# Pace et al. (2009) JBC 284:13285-13289, with non-ionizable pKs set
# to 7.4 (physiological pH)
K <- list(I = c( 7.4,   # Ala
                12.3,   # Arg
                 7.4,   # Asn
                 3.9,   # Asp
                 8.6,   # Cys
                 7.4,   # Gln
                 4.3,   # Glu
                 7.4,   # Gly
                 6.5,   # His
                 7.4,   # Ile
                 7.4,   # Leu
                10.4,   # Lys
                 7.4,   # Met
                 7.4,   # Phe
                 7.4,   # Pro
                 7.4,   # Ser
                 7.4,   # Thr
                 7.4,   # Trp
                 9.8,   # Tyr
                 7.4))  # Val
names(K$I) <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
                "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val")


# Given these biophysical indices, how similar are the amino acids? We have three-dimensions of measures here. Scatterplots can only display two dimensions ...
plot(Y$I, V$I, col="white", xlab = "hydrophobicity", ylab = "volume")
text(Y$I, V$I, names(Y$I))

plot(Y$I, K$I, col="white", xlab = "hydrophobicity", ylab = "pK")
text(Y$I, K$I, names(Y$I))

# ... but how do we plot 3D data? Plotting into a 3D cube is possible, but such
# plots are in general unintuitive and hard to interpret. One alternative is a
# so-called "ternary plot":

if (!require(ggtern)) {
  install.packages("ggtern")
  library(ggtern)
}
# Package information:
#  library(help = ggtern)       # basic information
#  browseVignettes("ggtern")    # available vignettes
#  data(package = "ggtern")     # available datasets



# collect into data frame, normalize to (0.05, 0.95)
myDat <- data.frame("phi" = 0.9*(((Y$I-min(Y$I))/(max(Y$I)-min(Y$I))))+0.05,
                    "vol" = 0.9*(((V$I-min(V$I))/(max(V$I)-min(V$I))))+0.05,
                    "pK"  = 0.9*(((K$I-min(K$I))/(max(K$I)-min(K$I))))+0.05,
                    stringsAsFactors = FALSE)
rownames(myDat) <- names(Y$I)

ggtern(data = myDat,
       aes(x = vol,
           y = phi,
           z = pK,
           label = rownames(myDat))) +
  geom_text()

# This results in a mapping of amino acids relative to each other that is
# similar to the Venn diagram you have seen in the notes.


# =    2  Mutation Data matrix  ================================================

# A mutation data matrix encodes all amino acid pairscores in a matrix.

# The Biostrings package contains the most common mutation data matrices.

if (!require(Biostrings, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biostrings")
  library(Biostrings)
}
# Package information:
#  library(help=Biostrings)       # basic information
#  browseVignettes("Biostrings")  # available vignettes
#  data(package = "Biostrings")   # available datasets

# Let's load the BLOSUM62 mutation data matrix from the package
data(BLOSUM62)

# ... and see what it contains. (You've seen this matrix before.)
BLOSUM62

# We can simply access values via the row/column names.
# Identical amino acids have high scores ...
BLOSUM62["H", "H"]   # Score for a pair of two histidines
BLOSUM62["S", "S"]   # Score for a pair of two serines

# Similar amino acids have low positive scores ...
BLOSUM62["L", "I"]   # Score for a leucine / lysine pair
BLOSUM62["F", "Y"]   # etc.

# Dissimilar amino acids have negative scores ...
BLOSUM62["L", "K"]   # Score for a leucine / lysine pair
BLOSUM62["Q", "P"]   # etc.


BLOSUM62["R", "W"]   # the matrix is symmetric!
BLOSUM62["W", "R"]


# =    3  Background score  ====================================================

# The mutation data matrix is designed to give high scores to homologous sequences, low scores to non-homologous sequences. What score on average should we expect for a random sequence?

# If we sample amino acid pairs at random, we will get a score that is the
# average of the individual pairscores in the matrix. Omitting the ambiguity
# codes and the gap character:

sum(BLOSUM62[1:20, 1:20])/400

# But that score could be higher for real sequences, for which the amino acid
# distribution is not random. For example membrane proteins have a large number
# of hydrophobic residues - an alignment of unrelated proteins might produce
# positive scores. And there are other proteins with biased amino acid
# compositions, in particular poteins that interact with multiple other
# proteins. Let's test how this impacts the background score by comparing a
# sequence with shuffled sequences. These have the same composition, but are
# obvioulsy not homologous. The data directory contains the FASTA file for the
# PDB ID 3FG7 - a villin headpiece structure with a large amount of
# low-complexity amino acid sequence ...

aa3FG7 <- readAAStringSet("./data/3FG7.fa")[[1]]

# ... and the FASTA file for the E. coli OmpG outer membrane porin (PDB: 2F1C)
# with an exceptionally high percentage of hydrophobic residues.

aa2F1C <- readAAStringSet("./data/2F1C.fa")[[1]]

# Here is a function that takes two sequences and
# returns their average pairscore.

averagePairScore <- function(a, b, MDM = BLOSUM62) {
  # Returns average pairscore of two sequences.
  # Parameters:
  #    a, b   chr   amino acid sequence string
  #    MDM          mutation data matrix. Default is BLOSUM62
  # Value:    num   average pairscore.
  a <- unlist(strsplit(a, ""))
  b <- unlist(strsplit(b, ""))
  v <- 0
  for (i in seq_along(a)) {
    v <- v + MDM[ a[i], b[i] ]
  }
  return(v / length(a))
}

orig3FG7 <- toString(aa3FG7)
orig2F1C <- toString(aa2F1C)
N <- 1000
scores3FG7 <- numeric(N)
scores2F1C <- numeric(N)
for (i in 1:N) {
  scores3FG7[i] <- averagePairScore(orig3FG7, toString(sample(aa3FG7)))
  scores2F1C[i] <- averagePairScore(orig2F1C, toString(sample(aa2F1C)))
}

# Plot the distributions
hist(scores3FG7,
     col="#5599EE33",
     breaks = seq(-1.5, 0, by=0.1),
     main = "Pairscores for randomly shuffled sequences",
     xlab = "Average pairscore from BLOSUM 62")
hist(scores2F1C,
     col="#55EE9933",
     breaks = seq(-1.5, 0, by=0.1),
     add = TRUE)
abline(v = sum(BLOSUM62[1:20, 1:20])/400, col = "firebrick", lwd = 2)
legend('topright',
       c("3FG7 (villin)", "2F1C (OmpG)"),
       fill = c("#5599EE33", "#55EE9933"), bty = 'n',
       inset = 0.1)

# This is an important result: even though we have shuffled significantly biased
# sequences, and the average scores trend above the average of the mutation data
# matrix, the average scores still remain comfortably below zero. This means
# that we can't (in general) improve a high-scoring alignment by simply
# extending it with randomly matched residues. We will only improve the score if
# the similarity of newly added residues is larger than what we expect to get by
# random chance!


# [END]
