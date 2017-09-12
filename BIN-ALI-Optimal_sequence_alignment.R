# BIN-ALI-Optimal_sequence_alignment.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-Optimal_sequence_alignment unit.
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

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...

# ==============================================================================

# = 1 Biostrings Pairwise Alignment

# Biostrings is one of the basic packages that the Bioconductor software
# landscape builds on. It stores sequences in "AAstring" objects and these are
# complex software structures that are designed to be able to handle
# genome-scale sequences. Biostrings functions - such as the alignment functions
# - expect their input to be Biostrings objects.
?AAString

AAString("ACDE")
s <- AAString("ACDE")
str(s)
# See: it's complicated. This is an "S4" object. Bioconductor uses these objects
# almost exclusively, but we will not be poking around in their internals. Just
# this: how do we get the sequence back out of an AAString object? The help page
# for XString - the parent "class" of AAStrings - mentions the  alternatives:

as.character(s)  # the base R version
toString(s)      # using the Biostrings function toString()

# While we need to remember to convert our sequences from the character vectors
# that we store in our database, to AAStrings that we can align, the alignment
# itself is really straightforward. The pairwiseAlignment() function was written
# to behave exactly like the functions you encountered on the EMBOSS server.

# First: make AAString objects ...
sel <- myDB$protein$name == "MBP1_SACCE"
aaMBP1_SACCE <- AAString(myDB$protein$sequence[sel])

sel <- myDB$protein$name == paste("MBP1_", biCode(YFO), sep = "")
aaMBP1_YFO <-   AAString(myDB$protein$sequence[sel])

?pairwiseAlignment

# ... and align.
# Global optimal alignment with end-gap penalties is default. (like EMBOSS needle)
ali1 <-  pairwiseAlignment(
  aaMBP1_SACCE,
  aaMBP1_YFO,
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10,
  gapExtension = 0.5)

str(ali1)  # Did you think the AAString object was complicated ?

# This is a Biostrings alignment object. But we can use Biostrings functions to
# tame it:
ali1
writePairwiseAlignments(ali1)   # That should look familiar

# And we can make the internal structure work for us  (@ is for classes as
# $ is for lists ...)
str(ali1@pattern)
ali1@pattern
ali1@pattern@range
ali1@pattern@indel
ali1@pattern@mismatch

# or work with "normal" R functions
# the alignment length
nchar(ali1@pattern)

# the number of identities
sum(s2c(as.character(ali1@pattern)) ==
      s2c(as.character(ali1@subject)))

# ... e.g. to calculate the percentage of identities
100 *
  sum(s2c(as.character(ali1@pattern)) ==
        s2c(as.character(ali1@subject))) /
  nchar(ali1@pattern)
# ... which should be the same as reported in the writePairwiseAlignments()
# output. Awkward to type? Then it calls for a function:
#
percentID <- function(al) {
  # returns the percent-identity of a Biostrings alignment object
  return(100 *
           sum(s2c(as.character(al@pattern)) ==
                 s2c(as.character(al@subject))) /
           nchar(al@pattern))
}

percentID(ali1)

# Compare with local optimal alignment (like EMBOSS Water)
ali2 <-  pairwiseAlignment(
  aaMBP1_SACCE,
  aaMBP1_YFO,
  type = "local",
  substitutionMatrix = "BLOSUM62",
  gapOpening = 50,
  gapExtension = 10)

writePairwiseAlignments(ali2)   # This has probably only aligned the N-terminal
# DNA binding domain - but that one has quite
# high sequence identity:
percentID(ali2)

# == TASK: ==

# Compare the two alignments. I have weighted the local alignment heavily
# towards an ungapped alignment by setting very high gap penalties. Try changing
# the gap penalties and see what happens: how does the number of indels change,
# how does the length of indels change...

# Fine. Please return to the Wiki to study BLAST alignment...


# ==============================================================================
#        PART FOUR: APSES Domain annotation by alignment
# ==============================================================================

# In this section we define the YFO APSES sequence by performing a global,
# optimal sequence alignment of the yeast domain with the full length protein
# sequence of the protein that was the most similar to the yeast APSES domain.
#

# I have annotated the yeast APSES domain as a proteinAnnotation in the
# database. To view the annotation, we can retrieve it via the proteinID and
# featureID. Here is the yeast protein ID:
myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"]

# ... assign it for convenience:
proID <- myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"]

# ... and if you look at the feature table, you can identify the feature ID
myDB$feature[ , c("ID", "name", "description")]
myDB$feature$ID[myDB$feature$name == "APSES fold"]

# ... assign it for convenience:
ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"]

# ... and with the two annotations we can pull the entry from the protein
# annotation table
myDB$proteinAnnotation[myDB$proteinAnnotation$protein.ID == proID &
                         myDB$proteinAnnotation$feature.ID == ftrID, ]

myDB$proteinAnnotation$ID[myDB$proteinAnnotation$protein.ID == proID &
                            myDB$proteinAnnotation$feature.ID == ftrID]

# ... assign it for convenience:
fanID <- myDB$proteinAnnotation$ID[myDB$proteinAnnotation$protein.ID == proID &
                                     myDB$proteinAnnotation$feature.ID == ftrID]

# The annotation record contains the start and end coordinates which we can use
# to define the APSES domain sequence with a substr() expression.
substr(myDB$protein$sequence[myDB$protein$ID == proID],
       myDB$proteinAnnotation$start[myDB$proteinAnnotation$ID == fanID],
       myDB$proteinAnnotation$end[myDB$proteinAnnotation$ID == fanID])

# Lots of code. But don't get lost. Let's recapitulate what we have done: we
# have selected from the sequence column of the protein table the sequence whose
# name is "MBP1_SACCE", and selected from the proteinAnnotation table the start
# and end coordinates of the annotation that joins an "APSES fold" feature with
# the sequence, and used the start and end coordinates to extract a substring.
# The expressions get lengthy, but it's not hard to wrap all of this into a
# function so that we only need to define name and feature.

dbGetFeatureSequence
dbGetFeatureSequence(myDB, "MBP1_SACCE", "APSES fold")


# Let's convert this to an AAstring and assign it:
aaMB1_SACCE_APSES <- AAString(dbGetFeatureSequence(myDB,
                                                   "MBP1_SACCE",
                                                   "APSES fold"))

# To align, we need the YFO sequence. Here is it's definition again, just
# in case ...

sel <- myDB$protein$name == paste("MBP1_", biCode(YFO), sep = "")
aaMBP1_YFO <- AAString(myDB$protein$sequence[sel])

# Now let's align these two sequences of very different length without end-gap
# penalties using the "overlap" type. "overlap" turns the
# end-gap penalties off and that is crucially important since
# the sequences have very different length.

aliApses <-  pairwiseAlignment(
  aaMB1_SACCE_APSES,
  aaMBP1_YFO,
  type = "overlap",
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10,
  gapExtension = 0.5)

# Inspect the result. The aligned sequences should be clearly
# homologous, and have (almost) no indels. The entire "pattern"
# sequence from QIYSAR ... to ... KPLFDF  should be matched
# with the "query". Is this correct?
writePairwiseAlignments(aliApses)

# If this is correct, you can extract the matched sequence from
# the alignment object. The syntax is a bit different from what
# you have seen before: this is an "S4 object", not a list. No
# worries: as.character() returns a normal string.
as.character(aliApses@subject)

# Now, what are the aligned start and end coordinates? You can read them from
# the output of writePairwiseAlignments(), or you can get them from the range of
# the match.

str(aliApses@subject@range)

# start is:
aliApses@subject@range@start

# ... and end is:
aliApses@subject@range@start + aliApses@subject@range@width - 1

# Since we have this section defined now, we can create a feature annotation
# right away and store it in myDB.  Copy the code-template below to your
# myCode.R file, edit it to replace the placeholder items with your data:
#
#  - The <PROTEIN ID> is to be replaced with the ID of MBP1_YFO
#  - The <FEATURE ID> is to be replaced with the ID of "APSES fold"
#  - <START> and <END> are to be replaced with the coordinates you got above
#
# Then execute the code and continue below the code template. If you make an
# error, there are instructions on how to recover, below.
#
# ===== begin code template: add a proteinAnnotation to the database =====

# == edit all placeholder items!
myProteinID <- "<PROTEIN ID>"
myFeatureID <- "<FEATURE ID>"
myStart <- <START>
  myEnd   <- <END>

  # == create the proteinAnnotation entry
  panRow <- data.frame(ID = dbAutoincrement(myDB$proteinAnnotation$ID),
                       protein.ID = myProteinID,
                       feature.ID = myFeatureID,
                       start = myStart,
                       end = myEnd,
                       stringsAsFactors = FALSE)
myDB$proteinAnnotation <- rbind(myDB$proteinAnnotation, panRow)

# == check that this was successful and has the right data
myDB$proteinAnnotation[nrow(myDB$proteinAnnotation), ]

# ===== end code template ===========================================

# ... continue here.
# I expect that a correct result would look something like
#          ID protein.ID feature.ID start end
# 63 my_fan_1   my_pro_1  ref_ftr_1     6 104

# If you made a mistake, simply overwrite the current version of myDB by loading
# your saved, good version:  load("myDB.01.RData") and correct your mistake

# If this is correct, save it
save(myDB, file = "myDB.02.RData")  # Note that it gets a new version number!

# Done with this part. Copy the sequence of the APSES domain of MBP1_<YFO> - you
# need it for the reverse BLAST search, and return to the course Wiki.



# = 1 Tasks




# [END]
