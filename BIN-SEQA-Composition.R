# BIN-SEQA-Composition.R
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-SEQA-Comparison unit
#
# Version: 1.1
#
# Date:    2017  11  -  2019  01
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
# Versions:
#           1.0    First live version 2017
#           0.1    First code copied from BCH441_A03_makeYFOlist.R
#
# TODO:
#
#
# == HOW TO WORK WITH LEARNING UNIT FILES ======================================
#
# DO NOT SIMPLY  source()  THESE FILES!
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
#  going on. That's not how it works ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                      Line
#TOC> ----------------------------------------------------------
#TOC>   1        Preparation                                  47
#TOC>   2        Aggregate properties                         68
#TOC>   3        Sequence Composition Enrichment             111
#TOC>   3.1        Barplot, and side-by-side barplot         134
#TOC>   3.2        Plotting ratios                           169
#TOC>   3.3        Plotting log ratios                       185
#TOC>   3.4        Sort by frequency                         200
#TOC>   3.5        Color by amino acid type                  215
#TOC> 
#TOC> ==========================================================================


# =    1  Preparation  =========================================================

if (! requireNamespace("seqinr", quietly = TRUE)) {
  install.packages("seqinr")
}
# Package information:
#  library(help = seqinr)       # basic information
#  browseVignettes("seqinr")    # available vignettes
#  data(package = "seqinr")     # available datasets

# Load a reference sequence to work with:

# If you have done the BIN-Storing_data unit:
   source("makeProteinDB.R")
   sel <- which(myDB$protein$name == sprintf("MBP1_%s", biCode(MYSPE)))
   mySeq <- myDB$protein$sequence[sel]

# If not, use the yeast Mbp1 sequence:
   mySeq <- dbSanitizeSequence(fromJSON("./data/MBP1_SACCE.json")$sequence)


# =    2  Aggregate properties  ================================================


# Let's try a simple function from seqinr: computing the pI of the sequence
?seqinr::computePI

# This takes as input a vector of upper-case AA codes

# We can use the function strsplit() to split the string
# into single characters

(s <- strsplit(mySeq, "")) # splitting on the empty spring
                           # splits into single characters
s <- unlist(s)             # strsplit() returns a list! Why?
                           # (But we don't need a list now...)

# Alternatively, seqinr provides
# the function s2c() to convert strings into
# character vectors (and c2s to convert them back).

seqinr::s2c(mySeq)


seqinr::computePI(s2c(mySeq))  # isoelectric point
seqinr::pmw(s2c(mySeq))        # molecular weight
seqinr::AAstat(s2c(mySeq))     # This also plots the distribution of
                               # values along the sequence

# A true Labor of Love has gone into the
# compilation of the "aaindex" data:

?aaindex
data(aaindex)  # "attach" the dataset - i.e. make it accessible as an
               # R object

length(aaindex)

# Here are all the index descriptions
for (i in 1:length(aaindex)) {
  cat(paste(i, ": ", aaindex[[i]]$D, "\n", sep=""))
}


# =    3  Sequence Composition Enrichment  =====================================


# Lets use one of the indices to calculate and plot amino-acid
# composition enrichment:
aaindex[[459]]$D

#
# Let's construct an enrichment plot to compare average frequencies
# with the amino acid counts in our sequence.

(refData <- aaindex[[459]]$I)                # reference frequencies in %
names(refData) <- seqinr::a(names(refData))  # change names to single-letter
                                             # code using seqinr's "a()" function
sum(refData)
refData        # ... in %


# tabulate the amino acid counts in mySeq
(obsData <- table(s2c(mySeq)))                # counts
(obsData <- 100 * (obsData / sum(obsData)))   # frequencies


# ==   3.1  Barplot, and side-by-side barplot  =================================

barplot(obsData, col = "#CCCCCC", cex.names = 0.7)
abline(h = 100/20, col="#BB0000")

barplot(refData, col = "#BB0000", cex.names = 0.7)
abline(h = 100/20, col="#555555")

# Ok: first problem - the values in obsData are in alphabetical order. But the
# values in refData are in alphabetical order of amino acid name: alanine,
# arginine, asparagine, aspartic acid ... A, R, N, D, E ... you will see this
# order a lot - one of the old biochemistry tropes in the field. So we need to
# re-order one of the vectors to match the other. That's easy though:
refData
(refData <- refData[names(obsData)])

barplot(refData, col = "#BB0000", cex.names = 0.7)
abline(h = 100/20, col="#555555")

# To compare the values, we want to see them in a barplot, side-by-side ...
barplot(rbind(obsData, refData),
        ylim = c(0, 12),
        beside = TRUE,
        col = c("#CCCCCC", "#BB0000"),
        cex.names = 0.7)
abline(h = 100/20, col="#00000044")

# ... and add a legend
legend (x = 1, y = 12,
        legend = c("mySeq", "Average composition"),
        fill = c("#CCCCCC", "#BB0000"),
        cex = 0.7,
        bty = "n")


# ==   3.2  Plotting ratios  ===================================================

# To better compare the values, we'll calculate ratios between
# obsData and refData

barplot(obsData / refData,
        col = "#CCCCCC",
        ylab = "Sequence / Average",
        cex.names = 0.7)
abline(h = 1, col="#BB0000")
abline(h = c(1/3, 3), lty = 2, col="#BB000055")

# ... but  ratios are not very good here, since the difference in height on the
# plot now depends on the order we compare in: ratios of 1/3 and 3 (dotted
# lines) are exactly the same fold-difference !

# ==   3.3  Plotting log ratios  ===============================================

# A better way to display this
# is to plot log(ratios).

barplot(log(obsData / refData),
        col = "#CCCCCC",
        ylab = "log(Sequence / Average)",
        cex.names = 0.7)
abline(h = log(1), col="#BB0000")
abline(h = log(c(1/3, 3)), lty = 2, col="#BB000055")

# Note how the three-fold difference lines are now the same distance from the
# line of equal ratio.

# ==   3.4  Sort by frequency  =================================================

barplot(sort(log(obsData / refData), decreasing = TRUE),
        ylim = c(-3.5,2),
        col = "#CCCCCC",
        ylab = "log(Sequence / Average)",
        cex.names = 0.7)
abline(h = log(1), col="#BB0000")
abline(h = log(c(1/3, 3)), lty = 2, col="#BB000055")

arrows(4, 1.8, 0, 1.8, length = 0.07)
text(5.5, 1.8, "Enriched", cex = 0.7)
arrows(20, 1.8, 24, 1.8, length = 0.07)
text(19.5, 1.8, "Depleted", pos = 2, cex = 0.7)

# ==   3.5  Color by amino acid type  ==========================================

# color the bars by type.
# define colors
AAcol <- character()
AAcol["A"] <- "#AABBAA"
AAcol["C"] <- "#FFEE77"
AAcol["D"] <- "#DD6600"
AAcol["E"] <- "#DD3300"
AAcol["F"] <- "#767D38"
AAcol["G"] <- "#BBBBCC"
AAcol["H"] <- "#A2A1FD"
AAcol["I"] <- "#70B6C6"
AAcol["K"] <- "#4563BB"
AAcol["L"] <- "#80C6B6"
AAcol["M"] <- "#AFCC34"
AAcol["N"] <- "#BB88CC"
AAcol["P"] <- "#7292B7"
AAcol["Q"] <- "#8866BB"
AAcol["R"] <- "#74A0FF"
AAcol["S"] <- "#9999CC"
AAcol["T"] <- "#99AADD"
AAcol["V"] <- "#9DB500"
AAcol["W"] <- "#76AD48"
AAcol["Y"] <- "#44CA97"

barplot(rep(1, 20), names.arg = names(AAcol), col = AAcol, cex.names = 0.5)

lR <- sort(log(obsData / refData), decreasing = TRUE)
barplot(lR,
        ylim = c(-3.5,2),
        col = AAcol[names(lR)],
        ylab = "log(Sequence / Average)",
        cex.names = 0.7)
abline(h = log(1), col="#00000055")
abline(h = log(c(1/3, 3)), lty = 2, col="#00000033")

arrows(4, 1.8, 0, 1.8, length = 0.07)
text(5.5, 1.8, "Enriched", cex = 0.7)
arrows(20, 1.8, 24, 1.8, length = 0.07)
text(19.5, 1.8, "Depleted", pos = 2, cex = 0.7)


# Task:
#   Interpret this plot. (Can you?) Which types of amino acids are enriched?
#   Depleted?




# [END]
