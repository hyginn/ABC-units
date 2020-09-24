# tocID <- "RPR-Biostrings.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Biostrings unit.
#
# Version:  1.2
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2020 Updates
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.0    2017 Revisions
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
#TOC>   Section  Title                                             Line
#TOC> -----------------------------------------------------------------
#TOC>   1        The Biostrings:: Package                            56
#TOC>   2        Getting Data into Biostrings:: Objects              88
#TOC>   3        Working with Biostrings:: Objects                  110
#TOC>   3.1        Properties                                       127
#TOC>   3.2        Subsetting                                       168
#TOC>   3.3        Operators                                        180
#TOC>   3.4        Transformations                                  187
#TOC>   4        Getting Data out of Biostrings:: Objects           194
#TOC>   5        More                                               203
#TOC>   5.1        Views                                            205
#TOC>   5.2        Iranges                                          219
#TOC>   5.3        StringSets                                       225
#TOC> 
#TOC> ==========================================================================


# This is a very brief introduction to the Biostrings:: package, other units will
# be using more of the Biostrings:: functions.


# =    1  The Biostrings:: Package  ============================================


# First, we install and load the Biostrings:: package from bioconductor (if we
# haven't done so already).

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
# Examine the package information:
library(help = Biostrings)       # basic information
browseVignettes("Biostrings")    # available vignettes
data(package = "Biostrings")     # available datasets


# At its core, Biostrings:: objects are "classes" of type XString (you can think
# of a "class" in R as a special kind of list), that can take on particular
# flavours for RNA, DNA or amino acid sequence information.

class(Biostrings::RNAString("AUG"))
class(Biostrings::DNAString("ATG"))
class(Biostrings::AAString("M"))

# An essential property of Biostrings:: objects is that they only allow letters
# from the applicable IUPAC alphabet:
Biostrings::RNAString("AUG")
Biostrings::DNAString("AUG")  # Error! No "U" in IUPAC DNA codes


# =    2  Getting Data into Biostrings:: Objects  ==============================


# Example: read FASTA. Extract sequence. Convert to DNAString object.
rawSeq <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")
rawSeq <- dbSanitizeSequence(rawSeq)
biosDNAseq <- Biostrings::DNAString(rawSeq) # converts the nucleotide sequence
                                            # into an object of class DNAstring

# Multi FASTA files can be read directly as a "XStringSet) ...
rawMFAfile <- "./data/S288C_YDL056W_MBP1_coding.fsa"
(biosDNASet <- Biostrings::readDNAStringSet(rawMFAfile))

# ... and if you subset one sequence from the set, you get an XString object
# back again.
(Xseq <- biosDNASet[[1]])

biosDNAseq == Xseq           # the comparison evaluates to TRUE ...
identical(biosDNAseq, Xseq)  # ... and indeed the objects are deemed identical.



# =    3  Working with Biostrings:: Objects  ===================================

# Biostrings:: is a highly engineered package that is tightly integrated into
# the Bioconductor world - unfortunately that brings with it a somewhat
# undesirable level of computational overhead and dependencies. Using the
# package as we normally do - i.e. calling required functions with their
# explicit package prefix is therefore not advisable. There are generics
# that won't be propery dispatched. If you only need a small number of
# functions for a very specific context, you will probably get away with
# Biostrings::<function>() - but even in the demonstration code of this script
# not everything works out of the box. We'll therefore load the library,
# but we'll (redundantly) use the prefix anyway so as to emphasize where
# the functions come from.

library(Biostrings)


# ==   3.1  Properties  ========================================================
str(rawSeq)
str(biosDNAseq)

length(rawSeq)       # ... is 1: one string only. To get the number of
                     # characters in a string, you need nchar().
length(biosDNAseq)   # but the length of a "Bstring" is the number of elements
nchar(rawSeq)
nchar(biosDNAseq)    # ... but nchar() works too.

(uL <- Biostrings::uniqueLetters(biosDNAseq))

# Count frequencies - with strings, you would strsplit() into a character
# vector and then use table(). biost
Biostrings::alphabetFrequency(biosDNAseq)

# letterFrequency() works with a defined alphabet - such as what uniqueLetters()
# returns.
Biostrings::letterFrequency(biosDNAseq, uL)
sum(Biostrings::letterFrequency(biosDNAseq, c("G", "C"))) /
  length(biosDNAseq) # GC contents

Biostrings::dinucleotideFrequency(biosDNAseq)
barplot(sort(Biostrings::dinucleotideFrequency(biosDNAseq)), cex.names = 0.5)

(triNuc <- Biostrings::trinucleotideFrequency(biosDNAseq))
barplot(sort(triNuc), col="#4499EE33")
triNuc[triNuc == max(triNuc)]
triNuc[triNuc == min(triNuc)]
max(triNuc) / min(triNuc)  # AAA is more than 13 times as frequent as CGT

# compare to a shuffled sequence:
(triNuc <- Biostrings::trinucleotideFrequency(sample(biosDNAseq)))
barplot(sort(triNuc), col="#EEEE4433", add = TRUE)
max(triNuc)
# Interpret this plot.
(triNuc <- Biostrings::trinucleotideFrequency(sample(biosDNAseq)))
barplot(sort(triNuc), col="#EEEE4433")
max(triNuc)


# ==   3.2  Subsetting  ========================================================

# Subsetting any XString object works as expected:
biosDNAseq[4:15]

# ... well - maybe not expected, because rawSeq[4:15] would not work.

# Alternatively to the "[" operator, use the subseq() function - especially for
# long sequences. This is far more efficient.
Biostrings::subseq(biosDNAseq, start = 1, end = 30)


# ==   3.3  Operators  =========================================================

# RNAstring() and DNAstring() objects compare U and T as equals!
  Biostrings::RNAString("AUGUCUAACCAAAUAUACUCAGCGAGAUAU") ==
  Biostrings::DNAString("ATGTCTAACCAAATATACTCAGCGAGATAT")


# ==   3.4  Transformations  ===================================================

biosDNAseq[4:15]
Biostrings::reverseComplement(biosDNAseq[4:15])
Biostrings::translate(biosDNAseq[4:15])


# =    4  Getting Data out of Biostrings:: Objects  ============================

# If you need a character object, use toString():

Biostrings::toString(biosDNAseq[4:15])

# saveRDS() and readRDS() works like on all other R objects.


# =    5  More  ================================================================

# ==   5.1  Views  =============================================================

# Biostring "Views" are objects that store multiple substrings of one
# Biostring object.

(myView <- Biostrings::Views(biosDNAseq,
                             start = c(1, 19, 37),
                             end = c(15, 30, 45)))

# Views are convenient to store feature annotations
names(myView) <- c("Feature-A", "Feature-B", "Feature-C")
cat(sprintf("\n%s\t(%d)\t%s", names(myView), width(myView), myView ))


# ==   5.2  Iranges  ===========================================================

# Biostrings:: Iranges are like Views with a common start point. These can be
# useful for feature annotations. Instead of start/end you store start/width.


# ==   5.3  StringSets  ========================================================

# Biostring "StringSets" store multiple sequences.
#
ompA <- Biostrings::AAString("MKKTAIAIAVALAGFATVAQA")
sample(ompA) # sample can work directly on a Biostring object to shuffle it

x <- Biostrings::toString(ompA)
for (i in 2:10) {
  x[i] <- Biostrings::toString(sample(ompA))
}
shuffledPeptideSet <- Biostrings::AAStringSet(x)
names(shuffledPeptideSet) <- c("ompA", paste("shuffle.", 1:9, sep=""))
shuffledPeptideSet

length(shuffledPeptideSet)
Biostrings::width(shuffledPeptideSet)
Biostrings::alphabetFrequency(shuffledPeptideSet)


# [END]
