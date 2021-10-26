# tocID <- "BIN-ALI-Optimal_sequence_alignment.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-Optimal_sequence_alignment unit.
#
# ==============================================================================
# Version:  1.7.1
#
# Date:     2017-09   -   2020-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.7.1  add jsonlite:: to fromjJSON() in code sample and ./myScripts/
#           1.7    2020 updates
#           1.6    Maintenance
#           1.5    Change from require() to requireNamespace(),
#                    use <package>::<function>() idiom throughout
#           1.4    Pull s2c() from seqinr package, rather then loading the
#                    entire library.
#           1.3    Updated confirmation task with correct logic
#           1.2    Added missing load of seqinr package
#           1.1    Update annotation file logic - it could already have been
#                    prepared in the BIN-FUNC-Annotation unit.
#           1.0.1  bugfix
#           1.0    First 2017 live version.
#           0.1    First code copied from 2016 material.
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
#TOC>   Section  Title                                                      Line
#TOC> --------------------------------------------------------------------------
#TOC>   1        Prepare                                                      58
#TOC>   2        Biostrings Pairwise Alignment                                75
#TOC>   2.1        Optimal global alignment                                   93
#TOC>   2.2        Optimal local alignment                                   156
#TOC>   3        APSES Domain annotation by alignment                        180
#TOC>   4        Update your database script                                 261
#TOC>   4.1        Preparing an annotation file ...                          267
#TOC>   4.1.1          If you HAVE NOT done the BIN-FUNC-Annotation unit     269
#TOC>   4.1.2          If you HAVE done the BIN-FUNC-Annotation unit         314
#TOC>   4.2        Execute and Validate                                      338
#TOC> 
#TOC> ==========================================================================


# =    1  Prepare  =============================================================

if (! requireNamespace("seqinr", quietly=TRUE)) {
  install.packages("seqinr")
}
# You can get package information with the following commands:
# library(help = seqinr)       # basic information
# browseVignettes("seqinr")    # available vignettes
# data(package = "seqinr")     # available datasets


# You need to recreate the protein database that you have constructed in the
# BIN-Storing_data unit.

source("./myScripts/makeProteinDB.R")


# =    2  Biostrings Pairwise Alignment  =======================================


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


# Biostrings stores sequences in "XString" objects. Once we have converted our
# target sequences to AAString objects, the alignment itself is straightforward.

# ==   2.1  Optimal global alignment  ==========================================

# The pairwiseAlignment() function was written to behave
# exactly like the functions you encountered on the EMBOSS server.

# First: make AAString objects ...
sel <- myDB$protein$name == "MBP1_SACCE"
aaMBP1_SACCE <- Biostrings::AAString(myDB$protein$sequence[sel])

sel <- myDB$protein$name == paste("MBP1_", biCode(MYSPE), sep = "")
aaMBP1_MYSPE <-   Biostrings::AAString(myDB$protein$sequence[sel])

?pairwiseAlignment
# ... and align.
# Global optimal alignment with end-gap penalties is default.
ali1 <-  Biostrings::pairwiseAlignment(
  aaMBP1_SACCE,
  aaMBP1_MYSPE,
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10,
  gapExtension = 0.5)

str(ali1)  # ... it's complicated

# This is a Biostrings alignment object. But we can use Biostrings functions to
# tame it:
ali1
Biostrings::writePairwiseAlignments(ali1)   # That should look familiar

# And we can make the internal structure work for us  (@ is for classes as
# $ is for lists ...)
str(ali1@pattern)
ali1@pattern
ali1@pattern@range
ali1@pattern@indel
ali1@pattern@mismatch

# or work with "normal" R functions
# the alignment length
nchar(as.character(ali1@pattern))

# the number of identities
sum(seqinr::s2c(as.character(ali1@pattern)) ==
    seqinr::s2c(as.character(ali1@subject)))

# ... e.g. to calculate the percentage of identities
100 *
  sum(seqinr::s2c(as.character(ali1@pattern)) ==
      seqinr::s2c(as.character(ali1@subject))) /
  nchar(as.character(ali1@pattern))
# ... which should be the same as reported in the writePairwiseAlignments()
# output. Awkward to type? Then it calls for a function:
#
percentID <- function(al) {
  # returns the percent-identity of a Biostrings alignment object
  return(100 *
         sum(seqinr::s2c(as.character(al@pattern)) ==
             seqinr::s2c(as.character(al@subject))) /
         nchar(as.character(al@pattern)))
}

percentID(ali1)

# ==   2.2  Optimal local alignment  ===========================================

# Compare with local optimal alignment (like EMBOSS Water)
ali2 <-  Biostrings::pairwiseAlignment(
  aaMBP1_SACCE,
  aaMBP1_MYSPE,
  type = "local",
  substitutionMatrix = "BLOSUM62",
  gapOpening = 50,
  gapExtension = 10)

Biostrings::writePairwiseAlignments(ali2)
# This has probably only aligned the N-terminal DNA binding domain - but that
# one has quite high sequence identity:
percentID(ali2)

# == TASK: ==

# Compare the two alignments. I have weighted the local alignment heavily
# towards an ungapped alignment by setting very high gap penalties. Try changing
# the gap penalties and see what happens: how does the number of indels change,
# how does the length of indels change...


# =    3  APSES Domain annotation by alignment  ================================

# In this section we define the MYSPE APSES sequence by performing a global,
# optimal sequence alignment of the yeast APSES domain with the full length
# protein sequence of the protein that was the most similar to the yeast APSES
# domain.
#

# I have annotated the yeast APSES domain as a feature in the
# database. To view the annotation, we can retrieve it via the proteinID and
# featureID. Here is the yeast protein ID:
(proID <- myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"])


# ... and if you look at the feature table, you can identify the feature ID
(ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"])

# ... and with the two annotations we can get the corresponding ID from the
# annotation table
(fanID <- myDB$annotation$ID[myDB$annotation$proteinID == proID &
                             myDB$annotation$featureID == ftrID])

myDB$annotation[myDB$annotation$ID == proID &
                myDB$annotation$ID == ftrID, ]

# The annotation record contains the start and end coordinates which we can use
# to define the APSES domain sequence with a substr() expression.

(start <- myDB$annotation$start[myDB$annotation$ID == fanID])
(end   <- myDB$annotation$end[myDB$annotation$ID == fanID])
(apses <- substr(myDB$protein$sequence[myDB$protein$ID == proID],
                 start,
                 end))

# Lots of code. But don't get lost. Let's recapitulate what we have done: we
# have selected from the sequence column of the protein table the sequence whose
# name is "MBP1_SACCE", and selected from the annotation table the start
# and end coordinates of the annotation that joins an "APSES fold" feature with
# the sequence, and used the start and end coordinates to extract a substring.

# Let's convert this to an AAstring and assign it:
aaMB1_SACCE_APSES <- Biostrings::AAString(apses)

# Now let's align these two sequences of very different length without end-gap
# penalties using the "overlap" type. "overlap" turns the
# end-gap penalties off and that is crucially important since
# the sequences have very different length.

aliApses <-  Biostrings::pairwiseAlignment(
  aaMB1_SACCE_APSES,
  aaMBP1_MYSPE,
  type = "overlap",
  substitutionMatrix = "BLOSUM62",
  gapOpening = 10,
  gapExtension = 0.5)

# Inspect the result. The aligned sequences should be clearly
# homologous, and have (almost) no indels. The entire "pattern"
# sequence from QIYSAR ... to ... KPLFDF  should be matched
# with the "query". Is this correct?
Biostrings::writePairwiseAlignments(aliApses)

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


# =    4  Update your database script  =========================================


# Since we have this feature defined now, we can create a feature annotation
# right away and store it in myDB.

# ==   4.1  Preparing an annotation file ...  ==================================
#
# ===   4.1.1  If you HAVE NOT done the BIN-FUNC-Annotation unit
#
#
#   You DON'T already have a file called "<MYSPE>-Annotations.json" in the
#   ./myScripts/ directory:
#
#   - Make a copy of the file "./data/refAnnotations.json" and put it in your
#     myScripts/ directory.
#
#   - Give it a name that is structured like "<MYSPE>-Annotations.json" - e.g.
#     if MYSPE is called "Crptycoccus neoformans", your file should be called
#     "CRYNE-Annotations.json" (and the "name" of your Mbp1 orthologue is
#     "MBP1_CRYNE").
#
#   - Open the file in the RStudio editor and delete all blocks for
#     the Mbp1 protein annotations except the first one.
#
#   - From that block, delete all lines except for the line that says:
#
# {"pName" : "MBP1_SACCE", "fName" : "APSES fold", "start" : "4", "end" : "102"},
#
#   - Then delete the comma at the end of the line (your file will just have
#     this one annotation).
#
#   - Edit that annotation: change MBP1_SACCE  to MBP1_<MYSPE> and change the
#     "start" and "end" features to the coordinates you just discovered for the
#     APSES domain in your sequence.
#
#   - Save the file in your myScripts/ directory
#
##   - Validate your file online at https://jsonlint.com/
#
#   - Update your "./myScripts/makeProteinDB.R" script to load your new
#     annotation when you recreate the database. Open the script in the
#     RStudio editor, and add the following command at the end:
#
#     myDB <- dbAddAnnotation(myDB,
#                 jsonlite::fromJSON("./myScripts/<MYSPE>-Annotations.json"))
#                                                 ^^^^^^^
#                                                edit this!
#   - save and close the file.
#
# Then SKIP the next section.
#
#
# ===   4.1.2  If you HAVE done the BIN-FUNC-Annotation unit    
#
#
#   You DO already have a file called "<MYSPE>-Annotations.json" in the
#   ./myScripts/ directory:
#
#   - Open the file in the RStudio editor.
#
#   - Below the last feature lines (but before the closing "]") add the
#     following feature line (without the "#")
#
# {"pName" : "MBP1_SACCE", "fName" : "APSES fold", "start" : "4", "end" : "102"}
#
#   - Edit that annotation: change MBP1_SACCE  to MBP1_<MYSPE> and change the
#     "start" and "end" features to the coordinates you just discovered for the
#     APSES domain in your sequence.
#
#   - Add a comma after the preceding feature line.
#
#   - Save your file.
#
#   - Validate your file online at https://jsonlint.com/
#
#
# ==   4.2  Execute and Validate  ==============================================
#
#   - source() your database creation script:
#
#  source("./myScripts/makeProteinDB.R")
#
#     This should run without errors or warnings. If it doesn't work and you
#     can't figure out quickly what's happening, ask on the mailing list for
#     help.
#
#   - Confirm
#     The following commands should retrieve the correct start and end
#     coordinates and sequence of the MBP1_MYSPE APSES domain:

sel <- which(myDB$protein$name == paste("MBP1_", biCode(MYSPE), sep = ""))

(proID <- myDB$protein$ID[sel])
(ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"])
(fanID <- myDB$annotation$ID[myDB$annotation$proteinID == proID &
                             myDB$annotation$featureID == ftrID])
(start <- myDB$annotation$start[myDB$annotation$ID == fanID])
(end   <- myDB$annotation$end[myDB$annotation$ID == fanID])
(apses <- substr(myDB$protein$sequence[myDB$protein$ID == proID],
                 start,
                 end))


# [END]
