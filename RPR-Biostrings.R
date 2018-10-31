# RPR-Biostrings.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Biostrings unit.
#
# Version:  1.0
#
# Date:     2017  10  20
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
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
#TOC>   Section  Title                                     Line
#TOC> ---------------------------------------------------------
#TOC>   1        The Biostrings Package                      52
#TOC>   2        Getting Data into Biostrings Objects        85
#TOC>   3        Working with Biostrings Objects            106
#TOC>   3.1      Properties                                 109
#TOC>   3.2      Subsetting                                 146
#TOC>   3.3      Operators                                  158
#TOC>   3.4      Transformations                            165
#TOC>   4        Getting Data out of Biostrings Objects     172
#TOC>   5        More                                       181
#TOC>   5.1      Views                                      183
#TOC>   5.2      Iranges                                    195
#TOC>   5.3      StringSets                                 201
#TOC>
#TOC> ==========================================================================


# This is a very brief introduction to the biostrings package, other units will
# be using more of the biostrings functions.


# =    1  The Biostrings Package  ==============================================


# First, we install and load the Biostrings package from bioconductor

if (! require(Biostrings, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("Biostrings")
  library(Biostrings)
}

# Examine the package information:
library(help = Biostrings)       # basic information
browseVignettes("Biostrings")    # available vignettes
data(package = "Biostrings")     # available datasets


# At its core, Biostrings objects are "classes" of type XString (you can think
# of a "class" in R as a special kind of list), that can take on particular
# flavours for RNA, DNA or amino acid sequence information.

class(RNAString("AUG"))
class(DNAString("ATG"))
class(AAString("M"))

# An essential property of Biostrings objects is that they only allow letters
# from the applicable IUPAC alphabet:
RNAString("AUG")
DNAString("AUG")  # Error! No "U" in IUPAC DNA codes


# =    2  Getting Data into Biostrings Objects  ================================


# Example: read FASTA. Extract sequence. Convert to DNAString object.
x <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")
x <- dbSanitizeSequence(x)
myDNAseq <- DNAString(x)   # takes the nucleotide sequence and converts into a
# object of class DNAstring

# Multi FASTA files can be read directly as a "XStringSet) ...
(myDNASet <- readDNAStringSet("./data/S288C_YDL056W_MBP1_coding.fsa"))

# ... and if you subset one sequence from the set, you get an XString object
# back again.
(Xseq <- myDNASet[[1]])

myDNAseq == Xseq           # the comparison evaluates to TRUE ...
identical(myDNAseq, Xseq)  # ... and indeed the objects are deemed identical.



# =    3  Working with Biostrings Objects  =====================================


# ==   3.1  Properties  ========================================================
str(myDNAseq)
length(myDNAseq)  # This gives you the _number of nucleotides_!
# By comparison ...
length(x)         # ... is 1: one string only. To get the number of
# characters in a string, you need nchar().
nchar(x)          # However ...
nchar(myDNAseq)   # ... also works.

uniqueLetters(myDNAseq)

# Count frequencies - with strings, you would strsplit() into a character
# vector and then use table(). biost
alphabetFrequency(myDNAseq)

# letterFrequency() works with a defined alphabet - such as what uniqueLetters()
# returns.
letterFrequency(myDNAseq, uniqueLetters(myDNAseq))

sum(letterFrequency(myDNAseq, c("G", "C"))) / length(myDNAseq) # GC contents

dinucleotideFrequency(myDNAseq)
barplot(sort(dinucleotideFrequency(myDNAseq)), cex.names = 0.5)

(triNuc <- trinucleotideFrequency(myDNAseq))
barplot(sort(triNuc), col="#4499EE33")
triNuc[triNuc == max(triNuc)]
triNuc[triNuc == min(triNuc)]
max(triNuc) / min(triNuc)  # AAA is more than 13 times as frequent as CGT

# compare to a shuffled sequence:
(triNuc <- trinucleotideFrequency(sample(myDNAseq)))
barplot(sort(triNuc), col="#EEEE4433", add = TRUE)

# Interpret this plot.


# ==   3.2  Subsetting  ========================================================

# Subsetting any XString object works as expected:
myDNAseq[4:15]

# ... well - maybe not expected, because x[4:15] would not work.

# Alternatively to the "[" operator, use the subseq() function - especially for
# long sequences. This is far more efficient.
subseq(myDNAseq, start = 1, end = 30)


# ==   3.3  Operators  =========================================================

# RNAstring() and DNAstring() objects compare U and T as equals!
RNAString("AUGUCUAACCAAAUAUACUCAGCGAGAUAU") ==
  DNAString("ATGTCTAACCAAATATACTCAGCGAGATAT")


# ==   3.4  Transformations  ===================================================

myDNAseq[4:15]
reverseComplement(myDNAseq[4:15])
translate(myDNAseq[4:15])


# =    4  Getting Data out of Biostrings Objects  ==============================

# If you need a character object, use toString():

toString(myDNAseq[4:15])

# save() and load() works like on all other R objects.


# =    5  More  ================================================================

# ==   5.1  Views  =============================================================

# Biostring "Views" are objects that store multiple substrings of one
# Biostring object.

(myView <- Views(myDNAseq, start = c(1, 19, 37), end = c(15, 30, 45)))

# Views are convenient to store feature annotations
names(myView) <- c("Feature-A", "Feature-B", "Feature-C")
cat(sprintf("\n%s\t(%d)\t%s", names(myView), width(myView), myView ))


# ==   5.2  Iranges  ===========================================================

# Biostrings Iranges are like Views with a common start point. These can be
# useful for feature annotations. Instead of start/end you store start/width.


# ==   5.3  StringSets  ========================================================

# Biostring "StringSets" store multiple sequences.
#
ompA <- AAString("MKKTAIAIAVALAGFATVAQA")
sample(ompA) # sample can work directly on a Biostring object to shuffle it

x[1] <- toString(ompA)
for (i in 2:10) {
  x[i] <- toString(sample(ompA))
}
shuffledPeptideSet <- AAStringSet(x)
names(shuffledPeptideSet) <- c("ompA", paste("shuffle.", 1:9, sep=""))
shuffledPeptideSet

length(shuffledPeptideSet)
width(shuffledPeptideSet)
alphabetFrequency(shuffledPeptideSet)


# [END]
