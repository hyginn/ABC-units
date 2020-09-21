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
# Version: 1.1
#
# Date:    2020-09-18
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.1    2020 Workflow changes
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
#TOC>   1        Preparations                      47
#TOC>   2        Suitable MYSPE Species            59
#TOC>   3        Adopt "MYSPE"                     83
#TOC> 
#TOC> ==========================================================================


# =    1  Preparations  ========================================================
#

# Execute the two conditionals below:
if (! file.exists("scripts/.myProfile.R")) {
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
# ./scripts/ABC-makeMYSPElist.R  In brief, data for genome-sequenced fungi
# was retrieved from https://fungi.ensembl.org; a search for homologues to
# yeast Mbp1 was performed with BLAST at the NCBI, and the data was merged.
# A representative organism at each genus-level was chosen from those hits
# that actual;ly have a homologue. Finally, a mapping table was constructed to
# asymmetrically retrieve unique species: a student number will retrieve
# a species, but (public) knowledge of the species cannot reconstruct the
# student number.

# Task: Study ./scripts/ABC-makeMYSPElist.R, it implements a typical workflow
#       of selecting and combining data from various data resources. Studying
#       it will give you a better sense of how such workflows can be
#       implemented in practice.


# =    3  Adopt "MYSPE"  =======================================================

# Execute:
( MYSPE <- getMYSPE(myStudentNumber) )

# If this produced an error, this session has not been properly set up. You
# may not yet have run  init()  and edited  .myProfile.R , or that file is not
# in your  myScripts/  folder. Fix this, and execute  source(".Rprofile") .

# If this produced NA, your Student Number may not be correct, or you are not
# in my class-list. Contact me.
# Otherwise, this should have printed a species name. Your unique species
# for this course.

biCode(MYSPE) # and what is it's "BiCode" ... ?

# Task: Note down the species name and its five letter BiCode on your Student
# Wiki user page. Use this species whenever this or future assignments refer
# to MYSPE. Whenever you start a session, it will automatically be loaded
# from  myScripts/.myProfile.R  and is available as  MYSPE .

# Here is some more information:
fungiDat <- read.csv("data/Species.csv")

# number of sequenced fungal genomes:
nrow(fungiDat)

# sequenced genomes of species:
sel <- MYSPE == gsub("^(\\S+\\s\\S+).*$", "\\1", fungiDat$Name)
( x <- fungiDat[sel, "Name"] )
Nspc <- length(x)  # save this for later ...

# sequenced genomes of genus:
sel <- gsub("\\s.*", "", MYSPE) == gsub("\\s.*", "", fungiDat$Name)
( x <- fungiDat[sel, "Name"] )
Ngen <- length(x) - Nspc

# order:
( x <- unique(fungiDat[sel, "Classification"]) )
Nord <- sum(fungiDat$Classification == x) - Ngen - Nspc
Nfng <- nrow(fungiDat) - Nord - Ngen - Nspc

# proportions
pCol <- c("#ed394e", "#ff9582", "#ffd5c4", "#f2f2f0")
pie(c(Nspc, Ngen, Nord, Nfng),
    labels = "",
    radius = 1,
    main = "MYSPE in genome-sequenced fungi",
    sub = MYSPE,
    col = pCol,
    clockwise = TRUE,
    init.angle = 90)
legend(x = 1.3, y = 0.8,     # position
       legend = c("Species", "Genus", "Order", "Fungi"),
       y.intersp = 1.5,      # line spacing for labels
       cex = 0.9,            # character size for labels
       bty = "n",            # "no" box around the legend
       pt.cex = 2,           # size of colour boxes
       pch = 15,
       col = pCol)

# [END]
