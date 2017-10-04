# BIN-SEQA-Comparison.R
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-SEQA-Comparison unit
#
# Version: 0.1
#
# Date:    2017  08  25
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 0.1    First code copied from BCH441_A03_makeYFOlist.R
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

# ==============================================================================
#        PART THREE: Sequence Analysis
# ==============================================================================


if (!require(seqinr, quietly=TRUE)) {
  install.packages("seqinr")
  library(seqinr)
}

help(package = seqinr) # shows the available functions

# Let's try a simple function
?computePI

# This takes as input a vector of upper-case AA codes
# Let's retrieve the MYSPE sequence from our datamodel
# (assuming it is the last one that was added):

db$protein[nrow(db$protein), "sequence"]

# We can use the function strsplit() to split the string
# into single characters

s <- db$protein[nrow(db$protein), "sequence"]
s <- strsplit(s, "") # splitting on the empty spring
# splits into single characters
s <- unlist(s)       # strsplit() returns a list! Why?
# (But we don't need a list now...)

# Alternatively, seqinr provides
# the function s2c() to convert strings into
# character vectors (and c2s to convert them back).

s <- s2c(db$protein[nrow(db$protein), "sequence"])
s

computePI(s)  # isoelectric point
pmw(s)        # molecular weight
AAstat(s)     # This also plots the distribution of
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


# Lets use one of the indices to calculate and plot amino-acid
# composition enrichment:
aaindex[[459]]

# === Sequence Composition Enrichment
#
# Let's construct an enrichment plot to compare one of the amino acid indices
# with the situation in our sequence.

refData <- aaindex[[459]]$I   # reference frequencies in %
names(refData) <- a(names(refData))  # change names to single-letter
# code using seqinr's "a()" function
refData


# tabulate our sequence of interest and normalize
obsData <- table(s)                # count occurrences
obsData = 100 * (obsData / sum(obsData))   # Normalize
obsData

len <- length(refData)

logRatio <- numeric() # create an empty vector

# loop over all elements of the reference, calculate log-ratios
# and store them in the vector
for (i in 1:len) {
  aa <- names(refData)[i] # get the name of that amino acid
  fObs <- obsData[aa]  # retrieve the frequency for that name
  fRef <- refData[aa]
  logRatio[aa] <- log(fObs / fRef) / log(2)  # remember log Ratio from
  # the lecture?
}

barplot(logRatio)

# Sort by frequency, descending
logRatio <- sort(logRatio, decreasing = TRUE)

barplot(logRatio)  # If you can't see all of the amino acid letters in the
# x-axis legend, make the plot wider by dragging the
# vertical pane-separator to the left



# label the y-axis
# (see text() for details)
label <- expression(paste(log[2],"( f(obs) / f(ref) )", sep = ""))

barplot(logRatio,
        main = paste("AA composition enrichment"),
        ylab = label,
        cex.names=0.9)



# color the bars by type.
# define colors
chargePlus  <- "#404580"
chargeMinus <- "#ab3853"
hydrophilic <- "#9986bf"
hydrophobic <- "#d5eeb1"
plain       <- "#f2f7f7"

# Assign the colors to the different amino acid names
barColors <- character(len)

for (i in 1:length(refData)) {
  AA <- names(logRatio[i])
  if (grepl("[HKR]",      AA)) {barColors[i] <- chargePlus }
  else if (grepl("[DE]",       AA)) {barColors[i] <- chargeMinus}
  else if (grepl("[NQST]",     AA)) {barColors[i] <- hydrophilic}
  else if (grepl("[FAMILYVW]", AA)) {barColors[i] <- hydrophobic}
  else                               barColors[i] <- plain
}

barplot(logRatio,
        main = paste("AA composition enrichment"),
        ylab = label,
        col = barColors,
        cex.names=0.9)


# draw a horizontal line at y = 0
abline(h=0)

# add a legend that indicates what the colours mean
legend (x = 1,
        y = -1,
        legend = c("charged (+)",
                   "charged (-)",
                   "hydrophilic",
                   "hydrophobic",
                   "plain"),
        bty = "n",
        fill = c(chargePlus,
                 chargeMinus,
                 hydrophilic,
                 hydrophobic,
                 plain)
)

# == TASK ==
# Interpret this plot. (Can you?)

#
#
# ==============================================================================








# [END]
