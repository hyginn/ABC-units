# tocID <- "BIN-ALI-Dotplot.R"
#
#
# ==============================================================================
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-Dotplot unit.
#
# Version:  0.2
#
# Date:     2019  01  07
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.2    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
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
#TOC>   Section  Title                  Line
#TOC> --------------------------------------
#TOC>   1        ___Section___            42
#TOC>   2        Tasks                   190
#TOC> 
#TOC> ==========================================================================


# =    1  ___Section___  =======================================================

if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly=TRUE)) {
  BiocManager::install("Biostrings")
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets

if (!requireNamespace("seqinr", quietly=TRUE)) {
  install.packages("seqinr")
}


# Let's load BLOSUM62
data(BLOSUM62, package = "Biostrings")

# Now let's craft code for a dotplot. That's surprisingly simple. We build a
# matrix that has as many rows as one sequence, as many columns as another. Then
# we go through every cell of the matrix and enter the pairscore we encounter
# for the amino acid pair whose position corresponds to the row and column
# index. Finally we visualize the matrix in a plot.
#

# First we fetch our sequences and split them into single characters.
sel <- myDB$protein$name == "MBP1_SACCE"
MBP1_SACCE <- seqinr::s2c(myDB$protein$sequence[sel])

sel <- myDB$protein$name == paste("MBP1_", biCode(MYSPE), sep = "")
MBP1_MYSPE <- seqinr::s2c(myDB$protein$sequence[sel])

# Check that we have two character vectors of the expected length.
str(MBP1_SACCE)
str(MBP1_MYSPE)

# How do we get the pairscore values? Consider: a single pair of amino acids can
# be obtained from sequence SACCE and MYSPE eg. from position 13 and 21 ...
MBP1_SACCE[13]
MBP1_MYSPE[21]

# ... using these as subsetting expressions, we can pull the pairscore
# from the MDM
BLOSUM62[MBP1_SACCE[13], MBP1_MYSPE[21]]

# First we build an empty matrix that will hold all pairscores ...
dotMat <- matrix(numeric(length(MBP1_SACCE) * length(MBP1_MYSPE)),
                 nrow = length(MBP1_SACCE), ncol = length(MBP1_MYSPE))

# ... then we loop over the sequences and store the scores in the matrix.
#
for (i in 1:length(MBP1_SACCE)) {
  for (j in 1:length(MBP1_MYSPE)) {
    dotMat[i, j] <- BLOSUM62[MBP1_SACCE[i], MBP1_MYSPE[j]]
  }
}

# Even though this is a large matrix, this does not take much time ...
# Let's have a look at a small block of the values:

dotMat[1:10, 1:10]

# Rows in this matrix correspond to an amino acid from MBP1_SACCE, columns in
# the matrix correspond to an amino acid from MBP1_MYSPE.

# To plot this, we use the image() function. Here, with default parameters.

image(dotMat)

# Be patient, this takes a few moments to render: more than 500,000 values.
# Nice.
# What do you expect?
# What would similar sequences look like?
# What do you see?

#You migh notice a thin line of yellow along the diagonal, moving approximately
# from bottom left to top right, fading in and out of existence. This is the
# signature of extended sequence similarity.

# Let's magnify this a bit by looking at only the first 200 amino acids ...
image(dotMat[1:200, 1:200])

# ... and, according to our normal writing convention, we would like the
# diagonal to run from top-left to bottom-right since we write from left to
# right and from top to bottom...
image(dotMat[1:200, 1:200], ylim = 1.0:0.0)

# ... and we would like the range of the x- and y- axis to correspond to the
# sequence position ...
image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1))

# ... and labels! Axis labels would be nice ...
image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1),
      xlab = "MBP1_MYSPE", ylab = "MBP1_SACCE" )

# ... and why don't we have axis-numbers on all four sides? Go, make that right
# too ...
len <- 200
image(x = 1:len, y = 1:len,  dotMat[1:len, 1:len], ylim=c(len,1),
      xlab = "MBP1_MYSPE", ylab = "MBP1_SACCE", axes = FALSE)
box()
axis(1, at = c(1, seq(10, len, by=10)))
axis(2, at = c(1, seq(10, len, by=10)))
axis(3, at = c(1, seq(10, len, by=10)))
axis(4, at = c(1, seq(10, len, by=10)))

# ... you get the idea, we can infinitely customize our plot. However a good way
# to do this is to develop a particular view for, say, a report or publication
# in a script and then put it into a function. I have put a function into the
# utilities file and called it dotPlot2(). Why not dotPlot() ... that's because
# there already is a dotplot function in the seqinr package:

seqinr::dotPlot(MBP1_SACCE, MBP1_MYSPE)                           # seqinr
dotPlot2(MBP1_SACCE, MBP1_MYSPE, xlab = "SACCE", ylab = "MYSPE")  # Our's

# Which one do you prefer? You can probably see the block patterns that arise
# from segments of repetitive, low complexity sequence. But you probably have to
# look very closely to discern the faint diagonals that correspond to similar
# sequence.


# Let's see if we can enhance the contrast between distributed noise and the
# actual alignment of conserved residues. We can filter the dot matrix with a
# pattern that enhances diagonally repeated values. Every value in the matrix
# will be replaced by a weighted average of its neighborhood. Here is  a
# diagonal-filter:

myFilter <- matrix(numeric(25), nrow = 5)
myFilter[1, ] <- c( 1, 0, 0, 0, 0)
myFilter[2, ] <- c( 0, 1, 0, 0, 0)
myFilter[3, ] <- c( 0, 0, 1, 0, 0)
myFilter[4, ] <- c( 0, 0, 0, 1, 0)
myFilter[5, ] <- c( 0, 0, 0, 0, 1)

# I have added the option to read such filters (or others that you could define on your own) as a parameter of the function.

dotPlot2(MBP1_SACCE, MBP1_MYSPE, xlab = "SACCE", ylab = "MYSPE", f = myFilter)

# I think the result shows quite nicely how the two sequences are globally
# related and where the regions of sequence similarity are. Play with this a bit
# ...  Can you come up with a better filter? If so, eMail us.




# =    2  Tasks  ===============================================================




# [END]
