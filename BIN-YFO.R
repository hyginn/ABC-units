# BIN-YFO.R
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-YFO unit
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
#        PART ONE: YFO Species
# ==============================================================================


# There are some "rabbitholes" that you are encouraged to follow to explore the code that went into generating the YFO species list. The minimal required result however is that you have picked an '''YFO''', and that its name got written into your personalized profile file.



# A detailed description of the process of compiling the YFO list of genome
# sequenced fungi with protein annotations and Mbp1 homologues is in the file
# ABC_makeYFOlist.R I encourage you to study it. Here we load the
# resulting vector of species name, then pick one of them in a random but
# reproducible way, determined by your student number.

load("data/YFOspecies.RData")  # load the species names
set.seed(myStudentNumber)      # seed the random number generator
YFO <- sample(YFOspecies, 1)   # pick a species at random
# write the result to your personalized profile data so we can use the result in
# other functions
cat(sprintf("YFO <- \"%s\"\n", YFO), file = ".myProfile.R", append = TRUE)

YFO         # so, which species is it ... ?
biCode(YFO) # and what is it's "BiCode" ... ?



# Note down the species name and its five letter label on your Student Wiki user page. '''Use this species whenever this or future assignments refer to YFO'''.



#
#
# ==============================================================================








# [END]
