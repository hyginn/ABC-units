# tocID <- "BIN-MYSPE.R"
#
# ---------------------------------------------------------------------------- #
#  PATIENCE  ...                                                               #
#    Do not yet work wih this code. Updates in progress. Thank you.            #
#    boris.steipe@utoronto.ca                                                  #
# ---------------------------------------------------------------------------- #
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-MYSPE unit
#
# Version: 1.0.1
#
# Date:    2017  09  21
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0.1  Move ABC-makeMYSPElist.R to ./scripts directory
# V 1.0    Final code, after rewriting BLAST parser and updating MYSPElist
# V 0.1    First code copied from BCH441_A03_makeMYSPElist.R
#
# TODO:
#
#
# == HOW TO WORK WITH LEARNING UNIT FILES ======================================
#
# DO NOT SIMPLY  source()  THESE FILES!
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
#  going on. That's not how it works ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                           Line
#TOC> -----------------------------------------------
#TOC>   1        Preparations                      39
#TOC>   2        Suitable MYSPE Species            51
#TOC>   3        Adopt "MYSPE"                     65
#TOC>
#TOC> ==========================================================================


# =    1  Preparations  ========================================================
#

# Execute the two conditionals below:
if (! file.exists(".myProfile.R")) {
  stop("PANIC: profile file does not exist. Fix problem or ask for help.")
}
if (! exists("myStudentNumber")) {
  stop("PANIC: profile data wasn't loaded. Fix problem or ask for help.")
}


# =    2  Suitable MYSPE Species  ==============================================


# In this unit we will select one species from a list of genome sequenced fungi
# and write it into your personalized profile file. This species will be called
# "MYSPE" (Your Favourite Organism) for other learning units and exercises.

# A detailed description of the process of compiling the list of genome
# sequenced fungi with protein annotations and Mbp1 homologues is in the file
# ./scripts/ABC-makeMYSPElist.R

# Task: Study ./scripts/ABC-makeMYSPElist.R, it implements a typical workflow
#       of selecting and combining data from  public-domain data resources.

# =    3  Adopt "MYSPE"  =======================================================


# In the code below, we load the resulting vector of species name, then pick one
# of them in a random but reproducible way, determined by your student number.

load("data/MYSPEspecies.RData")     # load the species names
set.seed(myStudentNumber)           # seed the random number generator
MYSPE <- sample(MYSPEspecies, 1)    # pick a species at random
set.seed(NULL)                      # reset the random number generator
# write the result to your personalized profile data so we can use the result in
# other functions
cat(sprintf("MYSPE <- \"%s\"\n", MYSPE), file = ".myProfile.R", append = TRUE)

MYSPE         # so, which species is it ... ?
biCode(MYSPE) # and what is it's "BiCode" ... ?

# Task: Note down the species name and its five letter label on your Student
# Wiki user page. Use this species whenever this or future assignments refer
# to MYSPE. In code, we will automatically load it from your .myProfile.R file.


# [END]
