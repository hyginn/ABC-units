# tocID <- "BIN-ALI-MSA.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-MSA unit.
#
# Version:  1.3
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.3    2020 updates
#           1.2    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
#           1.1    Added fetchMSAmotif()
#           1.0    Fully refactored and rewritten for 2017
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
#TOC>   Section  Title                                              Line
#TOC> ------------------------------------------------------------------
#TOC>   1        Preparations                                         55
#TOC>   2        Aligning full length MBP1 proteins                   97
#TOC>   2.1        Preparing Sequences                               108
#TOC>   2.2        Compute the MSA                                   133
#TOC>   3        Analyzing an MSA                                    154
#TOC>   4        Comparing MSAs                                      226
#TOC>   4.1        Importing an alignment to msa                     235
#TOC>   4.1.1          importing an .aln file                        244
#TOC>   4.1.2          Creating an MsaAAMultipleAlignment object     275
#TOC>   4.2        More alignments                                   326
#TOC>   4.3        Computing comparison metrics                      338
#TOC>   5        Profile-Profile alignments                          476
#TOC>   6        Sequence Logos                                      549
#TOC>   6.1        Subsetting an alignment by motif                  558
#TOC>   6.2        Plot a Sequence Logo                              607
#TOC> 
#TOC> ==========================================================================


# =    1  Preparations  ========================================================

# You need to reload you protein database, including changes that might
# have been made to the reference files. If you have worked with the
# prerequisite units, you should have a script named "./myScripts/makeProteinDB.R"
# that will create the myDB object with a protein and feature database.
# Ask for advice if not.
source("./myScripts/makeProteinDB.R")


if (! requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("Biostrings", quietly=TRUE)) {
  BiocManager::install("Biostrings")
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets

# Multiple sequence alignment algorithms are provided in
# the Bioconductor msa package.

if (! requireNamespace("msa", quietly=TRUE)) {
  BiocManager::install("msa")
}
# Package information:
#  library(help=msa)       # basic information
#  browseVignettes("msa")  # available vignettes
#  data(package = "msa")   # available datasets


# If an installation asks you if you want to update older packages, I recommend
# to always answer "a" for "all" unless you have an important reason not to. But
# if the installer asks you whether you want to compile from source, answer "n"
# for "no" unless you need new functionality in a particular bleeding-edge
# version of a package.

help(package = "msa")


# =    2  Aligning full length MBP1 proteins  ==================================

# In the Wiki part of this unit you have
#   - aligned full length MBP1 protein sequences at the EBI using T-Coffee
#   - saved the resulting alignment in CLUSTAL format
#     to the file "MBP1orthologuesTC.aln"

# In this section we will calculate an MSA of the same sequences using
# algorithms in the msa packages, and we will compare the two alignments.


# ==   2.1  Preparing Sequences  ===============================================

# We have used the biostrings' AAString() function before; for multiple
# alignments we need an AAStringSet(). AAStringSets are produced from vectors
# of sequence.

sel <- grep("MBP1", myDB$protein$name)
MBP1set <- Biostrings::AAStringSet(myDB$protein$sequence[sel])

# To help us make sense of the alignment we need to add the names for
# the sequences. Names for a seqSet object are held in the ranges slot...

MBP1set@ranges@NAMES <- myDB$protein$name[sel]

MBP1set

# You should have eleven sequences in this set, ask for advice if you don't.

# The little step of adding names is actually really very important. That's
# because the aligned sequences are meaningless strings of characters unless we
# can easily identify their biological relationships. Creating MSAs that are
# only identified by e.g. their RefSeq ids is a type of cargo-cult
# bioinformatics that we encounter a lot. The point of the alignment is not to
# create it, but to interpret it!

# ==   2.2  Compute the MSA  ===================================================

# The alignment itself is very simple. msa has msaMuscle() and msaClustalOmega()
# to produce alignments. (It also has msaClustalW() which is kind of
# embarrassing since that hasn't been the algorithm of first choice for
# decades. Don't use that one for real work.)


# Let's run an alignment with "Muscle"
(msaM <-  msa::msaMuscle( MBP1set, order = "aligned"))

# ... or to see the whole thing (cf. ?MsaAAMultipleAlignment ... print method):
msa::print(msaM, show=c("alignment", "complete"), showConsensus=FALSE)


# You see that the alignment object has sequence strings with hyphens as
# indel-characters. The names are printed to the console. And you also see that
# the order has not been preserved, but the most similar sequences are now
# adjacent to each other.


# =    3  Analyzing an MSA  ====================================================

# You probabaly realize that computing an MSA is not that hard. It's not
# entirely trivial to collect meaningful sequences via e.g. PSI-BLAST ... but
# then computing the alignment gives you a result quickly. But what does it
# mean? What information does the MSA contain?

# Let's have a first look at conserved vs. diverged regions of the MSA. msa
# provides the function msaConservationScore() which outputs a vector of scores.
# The scores are the sum of pairscores for the column: for example a perfectly
# conserved column of histidines would have the following score in our MSA
# of eleven sequences:
#   -  one (H, H) pair score is 8 in BLOSUM62;
#   -  there are (n^2 - n) / 2 pairs that can be formed between amino acids
#        in a column from n sequences;
#   -  therefore the column score is 8 * (11^2 - 11) / 2 == 440

# attach the BLOSUM62 package from the Biostrings package
data(BLOSUM62, package = "Biostrings")

msaMScores <- msa::msaConservationScore(msaM, substitutionMatrix = BLOSUM62)
plot(msaMScores, type = "l", col = "#205C5E", xlab = "Alignment Position")

# That plot shows the well-aligned regions (domains ?) of the sequence, but it
# could use some smoothing. Options for smoothing such plots include calculating
# averages in sliding windows ("moving average"), and lowess() smoothing. Here
# is a quick demo of a moving average smoothing, to illustrate the principle.

wRadius <- 15     # we take the mean of all values around a point +- wRadius
len <- length(msaMScores)
v <- msaMScores
for (i in (1 + wRadius):(len - wRadius)) {
  v[i] <- mean(msaMScores[(i - wRadius):(i + wRadius)]) # mean of values in
                                                        # window around i
}
points(v, col = "#FFFFFF", type = "l", lwd = 4.5)
points(v, col = "#3DAEB2", type = "l", lwd = 3)


# You can set a threshold and use rle() to define ranges of values that fall
# above and below the threshold, and thus approximate domain boundaries:
thrsh <- 30
(highScoringRanges <- rle(v > thrsh))
(idx <- c(1, cumsum(highScoringRanges$lengths)))
for (i in seq_along(highScoringRanges$lengths)) {
  if (highScoringRanges$values[i] == TRUE) { # If this range is above threshold,
    rect(idx[i], thrsh, idx[i+1], max(v),    # ... draw a rectangle
         col = "#205C5E33")                  # ... with a transparent color.
    cat(sprintf("Possible domain from %d to %d\n", idx[i], idx[i+1]))
  }
}

# Getting this right requires a bit of fiddling with the window radius and
# threshold (experiment with that a bit), but once we are satisfied, we can use
# the boundaries to print the MSA alignments separately for domains.
# Unfortunately the msa package provides no native way to extract blocks of the
# alignment for further processing, but your .utilities.R script contains a
# function to write alignment objects and sequence sets to .aln formatted
# output. Have a look:

writeALN

# Print out the aligned blocks
for (i in seq_along(highScoringRanges$lengths)) {
  if (highScoringRanges$values[i] == TRUE) { # If this range is above threshold,
    cat(sprintf("\n\nPossible domain from %d to %d\n", idx[i], idx[i+1]))
    writeALN(msaM, range = c(idx[i], idx[i+1]))
  }
}



# =    4  Comparing MSAs  ======================================================

# Let's compare the results of different alignment algorithms. We computed a
# vector of scores above, and we can compare that for different alignment
# algorithms. This is not trivial, so we'll need to look at that data in
# different ways and explore it. But first, let's get more alignments to compare
# with.


# ==   4.1  Importing an alignment to msa  =====================================

# We computed a T-Coffee alignment at the EBI. msa has no native import function
# so we need to improvise, and it's a bit of a pain to do - but a good
# illustration of strategies to convert data into any kind of object:
#   -  read an .aln file
#   -  adjust the sequence names
#   -  convert to msaAAMultipleAlignment object

# ===   4.1.1  importing an .aln file                   

# The seqinr package has a function to read CLUSTAL W formatted .aln files ...
if (! requireNamespace("seqinr", quietly=TRUE)) {
  install.packages("seqinr")
}
# Package information:
#  library(help=seqinr)       # basic information
#  browseVignettes("seqinr")  # available vignettes
#  data(package = "seqinr")   # available datasets

# read the T-coffee aligned file that you donwloaded from the EBI MSA tools
# (cf. http://steipe.biochemistry.utoronto.ca/abc/index.php/BIN-ALI-MSA)
tmp <- seqinr::read.alignment("msaT.aln", format = "clustal")

# read.alignment() returns a list. $seq is a list of strings, one for each
# complete alignment. However, they are converted to lower case.
x <- toupper(unlist(tmp$seq)) # get the sequences, uppercase

# $nam contains the names.
(names(x) <- tmp$nam)

#  Note that the names in the file are refseq IDs, we need to replace the
#  RefSeqIDs with the database names:
for (i in seq_along(x)) {
  # find the index of the RefSeqID
  id <- gsub("\\..*$", "", names(x)[i])     # fetch the name, drop the version
  sel <- which(myDB$protein$RefSeqID == id) # get the index
  names(x)[i] <- myDB$protein$name[sel]     # get the name
}

# ===   4.1.2  Creating an MsaAAMultipleAlignment object

# MsaAAMultipleAlignment objects are S4 objects that contain AAStringSet objects
# in their @unmasked slot, and a few additional items. Rather then build the
# object from scratch, we copy an existing object, and overwrite the data in its
# slots with what we need. Our goal is pragmatic, we want an object that msa's
# functions will accept as input.

# First: convert our named char vector into an AAstringSet
x <- Biostrings::AAStringSet(x)

# Then: create a new MsaAAMultipleAlignment S4 object. The msa package has
# defined what such an object should look like, with the SetClass() function. To
# create a new object, simply use the new() function, define the class that the
# object would have, and fill the slots with something that has the right type.
# How do we know the right type and the slot names? Easy! We just use
# str(<object>) to get the information.

str(msaM)

# There is a catch however in the way R makes such operations specific to
# the packages they need them: the function that creates the class is
# defined as a "generic", and when it is called, R looks in the package
# namespace for a more specific function with precise instructions what
# to do. However, we have not loaded the package namespace - we access all
# of the functions directly with the msa:: prefix. This method breaks down
# when generic functions are involved. I.e. - we could make it work, but
# the amount of code we need then is unreasonable. The straightforward
# way is to load the package. We can still use the prefix notation for
# its functions, just to emphasize where the function comes from. But since
# the namespace then exists, we ensure that generics are properly dispatched.

library(msa)  # load the msa package namespace

msaT <- new("MsaAAMultipleAlignment", # create new MsaAAMultipleAlignment object
            unmasked = x,             # "unmasked" slot takes an AAStringSet
            version = "T-Coffee",     # "version" slot takes a string
            params = list(),          # "params" takes a list(), we leave the
                                      #   list empty, but we could add the
                                      #   alignment parameters that we used at
                                      #   the EBI here.
            call = "imported from T-coffee alignment") # also a string

str(msaT)


msaT # Now we have fabricated an msaAAMultipleAlignment object, and we can
     # use the msa package functions on it

msaTScores <- msa::msaConservationScore(msaT, substitutionMatrix = BLOSUM62)

# ==   4.2  More alignments  ===================================================

# Next, we calculate alignments with msa's two other alignment options:
# CLUSTAL Omega
(msaO <- msa::msaClustalOmega( MBP1set, order = "aligned"))
msaOScores <- msa::msaConservationScore(msaO, substitutionMatrix = BLOSUM62)

# CLUSTAL W
(msaW <- msa::msaClustalW( MBP1set, order = "aligned"))
msaWScores <- msa::msaConservationScore(msaW, substitutionMatrix = BLOSUM62)


# ==   4.3  Computing comparison metrics  ======================================

# Ready to compare.

# ... sum of alignment scores of alignment divided by sum of alignment scores
#  of reference alignment (arbitrarily using CLUSTAL W as reference)

sRef <- sum(msaWScores)
sum(msaWScores) / sRef  # CLUSTAl W
sum(msaOScores) / sRef  # CLUSTAL O
sum(msaTScores) / sRef  # T-COFFEE
sum(msaMScores) / sRef  # MUSCLE

# ... mean alignment scores (higher is better)

mean(msaWScores)  # CLUSTAl W
mean(msaOScores)  # CLUSTAL O
mean(msaTScores)  # T-COFFEE
mean(msaMScores)  # MUSCLE

# total number of gaps (lower is better)
countGaps <- function(ali) {
  x <- paste0(as.character(ali), collapse = "")
  aa <- nchar(gsub("-", "", x))
  return(nchar(x) - aa)
}

countGaps(msaW)   # CLUSTAl W
countGaps(msaO)   # CLUSTAL O
countGaps(msaT)   # T-COFFEE
countGaps(msaM)   # MUSCLE

# number of indels in alignment (lower is less fragmented)
countIndels <- function(ali) {
  x <- paste0(as.character(ali), collapse = "")  # collapse into single string
  x <- unlist(strsplit(x, ""))                   # split into characters
  x <- x == "-"                                  # convert into boolean
  x <- rle(x)                                    # calculate rle
  # every run of TRUE is one indel event
  return(sum(x$values))
}

countIndels(msaW)   # CLUSTAl W
countIndels(msaO)   # CLUSTAL O
countIndels(msaT)   # T-COFFEE
countIndels(msaM)   # MUSCLE

# Let's look at the distribution of alignment scores:
boxplot(list(CLUSTAL.W = msaWScores,
             CLUSTAL.O = msaOScores,
             T.COFFEE  = msaTScores,
             MUSCLE    = msaMScores),
        cex.axis = 0.8,
        col = c("#7D556488", "#74628F88", "#5E78A188", "#3DAEB288"))

# CLUSTAL W and CLUSTAL O don't look all that different. T-Coffee tends to have
# a tighter distribution with less negative scores. Muscle has a slightly higher
# mean and generally higher scores.

# Boxplots are convenient, but don't give us much detail about the shape of the
# distribution. For that, we need histograms, or density plots.

plot(density(msaWScores),
     type = "l",
     col = "#7D5564",
     lwd = 1.5,
     ylim = c(0, (max(density(msaWScores)$y) * 1.3)),
     main = "Comparing MSA algorithms",
     xlab = "Alignment Score",
     ylab = "Density")
points(density(msaOScores),
       type = "l",
       lwd = 1.5,
       col = "#74628F")
points(density(msaTScores),
       type = "l",
       lwd = 1.5,
       col = "#5E78A1")
points(density(msaMScores),
       type = "l",
       lwd = 1.5,
       col = "#3DAEB2")
legend("topright",
       legend = c("MUSCLE", "T-COFFEE", "CLUSTAL O", "CLUSTAL W"),
       col = c("#3DAEB2", "#5E78A1", "#74628F", "#7D5564"),
       lwd = 2,
       cex = 0.7,
       bty = "n")

# The density plots confirm in more detail that CLUSTAL W misses some of the
# higher-scoring possibilities, that wherever CLUSTAL O is bad, it is quite bad,
# that T-COFFEE has fewer poorly scoring columns but misses some of the better
# scoring possibilities, and that MUSCLE appears to do best overall.

# Can we attribute these differences to sections of the alignment in which the
# algorithms did better or worse? Let's plot the scores cumulatively. The
# alignments have different lengths, so we plot the scores on the respective
# fraction of the alignement length.

plot(seq(0, 1, length.out = length(msaWScores)), # x- axis: fraction of length
     cumsum(msaWScores),
     type = "l",
     col = "#7D5564",
     lwd = 1.5,
     ylim = c(0, max(cumsum(msaMScores))),
     main = "Comparing MSA algorithms",
     xlab = "Alignment Position",
     ylab = "Cumulative Alignment Score")
points(seq(0, 1, length.out = length(msaOScores)), # x- axis: fraction of length
       cumsum(msaOScores),
       type = "l",
       lwd = 1.5,
       col = "#74628F")
points(seq(0, 1, length.out = length(msaTScores)), # x- axis: fraction of length
       cumsum(msaTScores),
       type = "l",
       lwd = 1.5,
       col = "#5E78A1")
points(seq(0, 1, length.out = length(msaMScores)), # x- axis: fraction of length
       cumsum(msaMScores),
       type = "l",
       lwd = 1.5,
       col = "#3DAEB2")
legend("bottomright",
       legend = c("MUSCLE", "T-COFFEE", "CLUSTAL O", "CLUSTAL W"),
       col = c("#3DAEB2", "#5E78A1", "#74628F", "#7D5564"),
       lwd = 2,
       cex = 0.7,
       bty = "n")

# Your alignment is going to be different from mine, due to the inclusion of
# MYSPE - but what I see is that MUSCLE gives the highest score overall, and
# achieves this with fewer indels than most, and the lowest number of gaps of
# all algorithms.

# To actually compare regions of alignments, we need to align alignments.


# =    5  Profile-Profile alignments  ==========================================


# Profile-profile alignments are the most powerful way to pick up distant
# relationships between sequence families. The can be used, for example to build
# a profile from structural superpositions of crystal structures, and then map a
# MSA alignment onto those features. Here we will use profile-profile comparison
# to compare two MSAs with each other, by aligning them. The algorithm is
# provided by the DECIPHER package.

if (! requireNamespace("DECIPHER", quietly=TRUE)) {
  BiocManager::install("DECIPHER")
}
# Package information:
#  library(help = DECIPHER)       # basic information
#  browseVignettes("DECIPHER")    # available vignettes
#  data(package = "DECIPHER")     # available datasets

# AlignProfiles() takes two AAStringSets as input. Let's compare the MUSCLE and
# CLUSTAL W alignments: we could do this directly ...
DECIPHER::AlignProfiles(msaW@unmasked, msaM@unmasked)

# But for ease of comparison, we'll reorder the sequences of the CLUSTAL W
# alignment into the same order as the MUSCLE alignment:
m <- as.character(msaM)
w <- as.character(msaW)[names(m)]

(ppa <- DECIPHER::AlignProfiles(AAStringSet(w), AAStringSet(m)))

# Conveniently, AlignProfiles() returns an AAStringSet, so we can use our
# writeALN function to show it. Here is an arbitrary block, from somewhere in
# the middle of the alignment:

writeALN(ppa, range = c(751, 810))

# If you look at this for a while, you can begin to figue out where the
# algorithms made different decisions about where to insert gaps, and how to
# move segments of sequence around. But matters become more clear if we
# post-process this profile-profile alignment. Let's replace all hyphens that
# the pp-alignment has inserted with blanks, and let's add a separator line down
# the middle between the two alignments.

x <- unlist(strsplit(as.character(ppa), ""))  # unlist all
dim(x) <- c(width(ppa)[1], length(ppa))       # form into matrix by columns
x <- t(x)                                     # transpose the matrix
(a1 <- 1:(nrow(x)/2))                         # rows of alignment 1
(a2 <- ((nrow(x)/2) + 1):nrow(x))             # rows of alignment 2
for (i in 1:ncol(x)) {
  if (all(x[a1, i] == "-")) { x[a1, i] <- " " }  # blank hyphens that shift
  if (all(x[a2, i] == "-")) { x[a2, i] <- " " }  # original alignment blocks
}

# collapse the matrix into strings
ppa2 <- character()
for (i in 1:nrow(x)) {
  ppa2[i] <- paste0(x[i, ], collapse = "")
  names(ppa2)[i] <- names(ppa)[i]
}

# add a separator line
x <- paste0(rep("-", width(ppa)[1]), collapse = "")
ppa2 <- c(ppa2[a1], x, ppa2[a2])

# inspect
writeALN(ppa2, range = c(800, 960))

# Again, go explore, and get a sense of what's going on. You may find that
# CLUSTAL W has a tendency to insert short gaps all over the alignment, whereas
# MUSCLE keeps indels in blocks. CLUSTAL's behaviour is exactly what I would
# expect from an algorithm that builds alignments incrementally from pairwise
# local alignments, without global refinement.


# =    6  Sequence Logos  ======================================================

# To visualize the information that we can get about structure and function with
# an MSA, we'll calculate a sequence logo of the Mbp1 recognition helix - the
# part of the structure that inserts into the major groove of the DNA and
# provides sequence specificity to the DNA binding of this transcription factor.
# Helix-B in Mbp1 with four residues upstream and downstream spans the sequence
#     43-ANFAKAKRTRILEKEVLKE-61

# ==   6.1  Subsetting an alignment by motif  ==================================

# Finding the location of such an substring in an alignment is not entirely
# trivial, because the alignment might have produced indels in that sequence.
# Our strategy can be:
#   -  fetch the sequence from the alignment
#   -  remove all hyphens
#   -  find the range where the target sequence matches
#   -  count how many characters in all there are in the aligned sequence, up
#        to the start and end of the match
#   -  these numbers define the range of the match in the alignment.

x <- as.character(msaM)["MBP1_SACCE"]
xAA <- gsub("-", "", x)

motif <- "ANFAKAKRTRILEKEVLKE"
(m <- regexpr(motif, xAA))  # matched in position 43, with a length of 19
(motifStart <- as.numeric(m))
(motifEnd <- attr(m, "match.length") + motifStart - 1)

# To count characters, we split the string into single characters ...
x <- unlist(strsplit(x, ""))

# ... convert this into a boolean, which is true if the character is not
# a hyphen ...
x <- x != "-"

# ... cast this into a numeric, which turns TRUE into 1 and FALSE into 0 ...
x <- as.numeric(x)

# ... and sum up the cumulative sum.
x <- cumsum(x)

# Now we can find where the 43'd and 61'st characters are located in the
# alignment string ...
(aliStart <- which(x == motifStart)[1])   # get the first hit if there are more
(aliEnd   <- which(x == motifEnd)[1])

# ... and subset the alignment

(motifAli <- subseq(msaM@unmasked, start = aliStart, end = aliEnd))

# Packaging this into a function is convenient to have, therefore I have added
# such a function to the .utilities.R script:  fetchMSAmotif(). Try it:

wing <- "HEKVQGGFGKYQGTWV" # the MBP1_SACCE APSES "wing" sequence
writeALN(fetchMSAmotif(msaM, wing))


# ==   6.2  Plot a Sequence Logo  ==============================================

# The Bioconductor seqLogo:: packager handles only DNA sequences, but there are
# now several good options to plot sequence logos in R, these include dagLogo,
# DiffLogo, Logolas, and motifStack. For our example we will use ggseqlogo
# written by by Omar Waghi, a former UofT BCB student who is now at the EBI.

if (! requireNamespace("ggseqlogo", quietly=TRUE)) {
  install.packages(("ggseqlogo"))
}
# Package information:
#  library(help=ggseqlogo)       # basic information
#  browseVignettes("ggseqlogo")  # available vignettes
#  data(package = "ggseqlogo")   # available datasets

ggseqlogo::ggseqlogo(as.character(motifAli))





# [END]
