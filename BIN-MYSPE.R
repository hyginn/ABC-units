# tocID <- "BIN-MYSPE.R"
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-MYSPE unit
#
#
# Version: 1.4.1
#
# Date:    2017-09  -  2022-10
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.4.1  Remove Wiki task - record MYSPE in the "About me" document instead
# V 1.4    Add troubleshooting hints via errText[[...]]
# V 1.3    2021 update of MYSPE mechanics; fix a bug no one had complained about
# V 1.2    Reorganized proportional plot section into a "further reading"
#          section, added nested-box, and sankey plot visualization of
#          proportions. Introduced plotly.
# V 1.1    2020 Workflow changes
# V 1.0.1  Move ABC-makeMYSPElist.R to ./scripts directory
# V 1.0    Final code, after rewriting BLAST parser and updating MYSPElist
# V 0.1    First code copied from BCH441_A03_makeMYSPElist.R
#
# TODO:    Sample solution for sankey plot function.
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
#TOC>   Section  Title                                             Line
#TOC> -----------------------------------------------------------------
#TOC>   1        PREPARATIONS                                        52
#TOC>   2        SUITABLE MYSPE SPECIES                              65
#TOC>   3        ADOPT "MYSPE"                                       89
#TOC>   4        FURTHER READING: PLOTTING PROPORTIONS              130
#TOC>   4.1        Percentages                                      148
#TOC>   4.2        Visualizing proportions: Pie chart               167
#TOC>   4.3        Visualizing proportions: Nested squares          245
#TOC>   4.4        Visualizing proportions: Sankey diagrams         282
#TOC>
#TOC> ==========================================================================


# =    1  PREPARATIONS  ========================================================
#

# Execute the two conditionals below:
if (! file.exists("./myScripts/.myProfile.R")) {
  stop(errText[["noProfileFile"]])     # message defined in .Rprofile
}

if (! exists("myStudentNumber")) {
  stop(errText[["noStudentNumber"]])   # message defined in .Rprofile
}


# =    2  SUITABLE MYSPE SPECIES  ==============================================


# In this unit we will select one species from a list of genome sequenced fungi
# and write it into your personalized profile file. This species will be called
# "MYSPE" (My Species) for other learning units and exercises.

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


# =    3  ADOPT "MYSPE"  =======================================================

# Execute:
( MYSPE <- getMYSPE(myStudentNumber) )

# If this produced an error, this session has not been properly set up. You
# may not yet have run  init()  and edited  .myProfile.R , or that file is not
# in your  myScripts/  folder. Fix this, and execute:
#
#    source(".Rprofile") .

# If this produced NA, your Student Number may not be correct, or you are not in
# my class-list. Contact me. Otherwise, this should have printed a species name,
# and the taxonomy ID of its genome-sequenced strain. This is your unique
# speciesfor this course. Note it in your journal ...

biCode(MYSPE) # and also note it's "BiCode" ...
( myTaxID <- names(MYSPE) )  # and its taxID


# Task:
# =====

#   Note down the species name and its five letter BiCode in the "About me"
#   document in your shared course folder. Use this species whenever this or
#   future assignments refer to MYSPE. Whenever you start a session, it will
#   automatically be loaded from  myScripts/.myProfile.R  and is available as
#   MYSPE .

# Here is some more information about MYSPE, taken from the table of genome-
# sequenced fungi that is in your ./data folder.
fungiDat <- read.csv("data/Species.csv")
iMs <- which(fungiDat$Taxon.ID == myTaxID)

( myOr <- fungiDat$Classification[iMs] )  # Taxonomic order
( myGn <- gsub("\\s.*", "", MYSPE))       # Taxonomic genus
( mySt <- fungiDat$Name[iMs] )            # Taxonomic strain

# That's all.


# =    4  FURTHER READING: PLOTTING PROPORTIONS  ===============================

# The material below is an exploration of data-preparation and plotting
# techniques; you can treat this as additional practice and further reading and
# I expect that some of the code and plotting examples may be useful in a
# different context.

# A frequent task is to visualize the proportion of elements with given
# categories in a sample. For example, we might ask what the proportion of the
# different orders of fungi is the order of MYSPE? Let's first collect the
# numbers.

( nFungi <- nrow(fungiDat) )                            # sequenced fungi
( nOrder <- sum(grepl(myOr, fungiDat$Classification)) ) # same order as MYSPE
( nGenus <- sum(grepl(myGn, fungiDat$Name)) )           # same genus as MYSPE
( nSpecies <- sum(grepl(MYSPE, fungiDat$Name)) )        # same species as MYSPE


# ==   4.1  Percentages  =======================================================

# The zeroth-order approach to visualization is simply to print percentages:

cat(sprintf("\n%s comprise %5.2f%% of fungi.",
        myOr,
        (nOrder * 100) / nFungi))

# ... or, adding the actual numbers:

cat(sprintf("\n%s comprise %5.2f%% of fungi (%d of %d).",
            myOr,
            (nOrder * 100) / nFungi,
            nOrder,
            nFungi))

# But that's hard to visualize for most of us, and anyway, we don't know how
# that relates to other orders.

# ==   4.2  Visualizing proportions: Pie chart  ================================

# Often, we will use a pie chart instead. Pie charts are rather informal types
# of plots, not well suited for analysis. But easy to do:

# Define four colors to identify the four categories
pCol <- c("#ed394e", "#ff9582", "#ffd5c4", "#f2f2f0")

oPar <- par(mar = c(0.1, 0.1, 2.5, 0.1))   # set margins to ~ 0
                                           # and remember the
                                           # previous setting

pie(c(nSpecies,                            # subtract numbers since these
      nGenus - nSpecies,                   # categories are mutually contained
      nOrder - nGenus - nSpecies,          # in each other
      nFungi - nOrder - nGenus - nSpecies),
      labels = "",
      radius = 0.9,
      main = "MYSPE in genome-sequenced fungi",
      lty = 0,                             # turn borders for wedges off
      col = pCol,
      clockwise = TRUE,
      init.angle = 90)

title(main=MYSPE, line=0, cex.main=0.7)    # add a title to the plot

legend(x = 0.95, y = 0.8,    # place at legend here
       legend = c("Species", "Genus", "Order", "Fungi"),
       y.intersp = 2,                      # line spacing for labels
       cex = 0.8,                          # character size for labels
       bty = "n",                          # "no" box around the legend
       pt.cex = 2,                         # size of colour boxes
       pch = 15,                           # a filled square
       col = pCol)

par(oPar)                                  # reset graphics state

# Unless MYSPE is one of the frequently sequenced species, there will only be a
# very thin wedge visible. Pie charts are not well suited to visualize small
# proportions.

# It is a little more useful if we have non-nested proportions - like the
# number of species in the same order overall:

myTbl <- sort(table(fungiDat$Classification), decreasing = TRUE)
head(myTbl)

# pie() does a reasonable job out of the box to interpret table() data:
pie(myTbl)

# ... we can improve this quickly with a bit of tweaking:

N <- length(myTbl)
sel <- myOr == names(myTbl) # TRUE for the MYSPE order, FALSE elsewhere

myCol <- rep(pCol[4], N)       # N elements of pCol[1]
myCol[sel] <- pCol[1]          # replace this one color

myLbl <- rep("", N)            # N labels of ""
myLbl[sel] <- myOr             # replace this one label with the MYSPE order


oPar <- par(mar = c(0.1, 0.1, 2.5, 0.1))   # set margins to ~ 0

pie(myTbl,
    labels = myLbl,
    radius = 0.9,
    main = "MYSPE order",
    border = "#DDDDDD",
    col = myCol,
    clockwise = TRUE,
    init.angle = 90)

par(oPar)                                  # reset graphics state

# But the overall problem remains.


# ==   4.3  Visualizing proportions: Nested squares  ===========================

# A simple alternative is to draw such proportions as nested squares:

x <- sqrt(nFungi)

# set margins to ~ 0 and type to square
oPar <- par(mar = c(0.1, 0.1, 0.1, 0.1), pty = "s")

# empty, square plot
plot(c(0, x), c(0, x), xlim = c(0, x), ylim = c(0, x),
     type="n", axes=FALSE, xlab="", ylab="")

# basic square for all genomes
rect(0, 0, x,              x,              col = pCol[4])

# grid
u <- 0:floor(x)
N <- length(u)
segments(rep(0, N), u, rep(x, N), u, col = "#0000FF18")
segments(u, rep(0, N), u, rep(x, N), col = "#0000FF18")
# each square on this grid is one genome

# colored squares
rect(0, 0, sqrt(nOrder),   sqrt(nOrder),   col = pCol[3])
rect(0, 0, sqrt(nGenus),   sqrt(nGenus),   col = pCol[2])
rect(0, 0, sqrt(nSpecies), sqrt(nSpecies), col = pCol[1])

# labels
text(x/2, x/2,      "Fungi")
text(x * 0.08, x * 0.11, myOr,   pos = 4, cex = 0.9)
text(x * 0.08, x * 0.06, myGn,   pos = 4, cex = 0.8)
text(x * 0.08, x * 0.02, MYSPE, pos = 4, cex = 0.7)

par(oPar)                                  # reset graphics state


# ==   4.4  Visualizing proportions: Sankey diagrams  ==========================

# Sankey diagrams are an excellent way to visualize complicated nested
# proportions and their changes (see here for example:
# https://www.r-graph-gallery.com/sankey-diagram.html). Here is a very simple
# example with the MYSPE proportions, as an illustration of the plotting
# principle.

if (! requireNamespace("plotly")) {
  install.packages("plotly")
}
# Package information:
#  library(help   = plotly)     # basic information
#  browseVignettes("plotly")    # available vignettes
#  data(package  = "plotly")    # available datasets

# Here, we use the plotly package that wraps a very well developed javascript
# library with many options for interactive plots. I am producing this plot
# hard-coded for the sample organism "Sporothrix schenkii"; you would need
# to change the code to adapt it to your own MYSPE - or even build a function
# for this. Do try this if you have a bit of coding experience, sankey diagrams
# are a good way to show hierarchical data relations - and if you get this
# working for your own organism you can be proud that you have understood
# how preparing the data works.


myNodes <- list(label = c("Fungi (1014)",              # 0 <- node ID
                          "Ophiostomatales (6)",       # 1
                          "Other...",                  # 2
                          "Sporothrix (4)",            # 3
                          "Other...",                  # 4
                          "Sporothrix schenckii (2)",  # 5
                          "Other..."                   # 6
                          ),
                x = c(0.1, 0.4, 0.4, 0.7, 0.7, 1.0, 1.0),
                y = c(0.3, 0.1, 0.7, 0.2, 0.7, 0.3, 0.7),
                color = c("#f2f2f0", #
                          "#ffd5c4",
                          "#CCCCCC",
                          "#ff9582",
                          "#CCCCCC",
                          "#ed394e",
                          "#CCCCCC"
                          ),
                pad = 15,
                thickness = 20,
                line = list(color = "black",
                            width = 0.5))

myLinks <- list(source = c(0, 0, 1, 1, 3, 3),   # i.e. there is a link of
                target = c(1, 2, 3, 4, 5, 6),   # weight 6 between node 0
                value =  c(6, 18, 4, 2, 2, 2))  # and node 1

# Setting up the actual plot ...
fig  <-  plotly::plot_ly(type = "sankey",
                         arrangement = "snap",
                         orientation = "h",
                         node = myNodes,
                         link = myLinks)

# Adding and adjusting a few layout parameters
fig <- plotly::layout(fig,
              title = "Fungi Genomes - Classification",
              font = list(size = 10))

fig     # plot the diagram

# Note that the plot appears in the Viewer window, not the Plot window, and that
# it is interactive: you can hover over nodes and links, and drag the nodes
# around.

# [END]
