# BIN-FUNC-Domain_annotation.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-FUNC-Domain_annotation unit.
#
# Version:  1.0
#
# Date:     2017  11 13
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    Live version 2017
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
#TOC>   Section  Title                                                               Line
#TOC> -----------------------------------------------------------------------------------
#TOC>   1        Update your database script                                           41
#TOC>   1.1      Preparing an annotation file ...                                      47
#TOC>   1.1.1    If you HAVE NOT done the BIN-ALI-Optimal_sequence_alignment unit      49
#TOC>   1.1.2    If you HAVE done the BIN-ALI-Optimal_sequence_alignment               93
#TOC>   1.2      Execute and Validate                                                 119
#TOC>   2        Plot Annotations                                                     144
#TOC> 
#TOC> ==========================================================================


# =    1  Update your database script  =========================================


# Since you have recorded domain features at the SMART database, we can store
# the feature annotations in myDB.

# ==   1.1  Preparing an annotation file ...  ==================================
#
# ===  1.1.1  If you HAVE NOT done the BIN-ALI-Optimal_sequence_alignment unit
#
#
#   You DON'T already have a file called "<MYSPE>-Annotations.json" in the
#   ./data/ directory:
#
#   - Make a copy of the file "./data/refAnnotations.json" and put it in your
#     project directory.
#
#   - Give it a name that is structured like "<MYSPE>-Annotations.json" - e.g.
#     if MYSPE is called "Crptycoccus neoformans", your file should be called
#     "CRYNE-Annotations.json" (and the "name" of your Mbp1 orthologue is
#     "MBP1_CRYNE").
#
#   - Open the file in the RStudio editor and delete all blocks for
#     the Mbp1 protein annotations except the first one.
#
#   - From that block, delete all lines that have annotations you did not
#     find in SMART for MBP1_MYSPE.
#
#   - Make enough copies of the "Ankyrin fold" and "low complexity" region
#     lines to have a line for each feature you found.
#
#   - Then delete the comma at the end of the last line.
#
#   - Edit the annotations: change MBP1_SACCE  to MBP1_<MYSPE> everywhere
#     and change the "start" and "end" features to the coordinates you
#     recorded in the SMART database.
#
#   - Save your file.
#
#   - Validate your file online at https://jsonlint.com/
#
#   - Update your "makeProteinDB.R" script to load your new
#     annotation when you recreate the database. Open the script in the
#     RStudio editor, and add the following command at the end:
#
#     myDB <- dbAddAnnotation(myDB, fromJSON("<MYSPE>-Annotations.json"))
#
#   - save and close the file.
#
# Then SKIP the next section.
#
#
# ===  1.1.2  If you HAVE done the BIN-ALI-Optimal_sequence_alignment         
#
#
#   You DO already have a file called "<MYSPE>-Annotations.json" in the
#   ./data/ directory:
#
#   - Open the file in the RStudio editor.
#
#   - Make as many copies of the "APSES fold" line as you have found
#     features in SMART.
#
#   - Add a comma after every line except for the last one
#
#   - Edit the annotations but include only features that are in the
#     myDB$feature table. Check which features are in the databse by executing
#
#        myDB$feature$name
#
#   - Update the "start" and "end" coordinates for each feature to the
#     values you found.
#
#   - Save your file.
#
#   - Validate your file online at https://jsonlint.com/
#
#
# ==   1.2  Execute and Validate  ==============================================
#
#   - source() your database creation script:
#
#     source("makeProteinDB.R")
#
#     This should run without errors or warnings. If it doesn't work and you
#     can't figure out quickly what's happening, ask on the mailing list for
#     help.
#
#   - Confirm
#     The following commands should retrieve all of the features that have been
#     annotated for MBP1_MYSPE

sel <- myDB$protein$name == paste("MBP1_", biCode(MYSPE), sep = "")

(proID  <- myDB$protein$ID[sel])
(fanIDs <- myDB$annotation$ID[myDB$annotation$proteinID == proID])
(ftrIDs <- unique(myDB$annotation$featureID[fanIDs]))
myDB$feature$name[ftrIDs] # This should list ALL of your annotated features
                          # (once). If not, consider what could have gone wrong
                          # and ask on the list if you have difficulties fixing
                          # it.


# =    2  Plot Annotations  ====================================================

# In this section we will plot domain annotations as colored rectangles on a
# sequence, as an example for using the R plotting system for generic, data
# driven images.

# We need a small utility function that draws the annotation boxes on a
# representation of sequence. It should accept the start and end coordinates,
# the y value where it should be plotted and the color of the box, and plot a
# rectangle using R's rect() function.

drawBox <- function(xStart, xEnd, y, myCol) {
  # Draw a box from xStart to xEnd at y, filled with colour myCol
  delta <- 0.1
  rect(xStart, (y - delta), xEnd, (y + delta),
       border = "black", col = myCol)
}

# test this:
plot(c(-1.5, 1.5), c(0, 0), type = "l")
drawBox(-1, 1, 0.0, "peachpuff")

# Next, we define a function to plot annotations for one protein: the name of
# the protein, a horizontal grey line for its length, and all of its features.

plotProtein <- function(DB, name, y) {
  # DB: protein database
  # name: the name of the protein in the database.
  # y: height where to draw the plot
  #
  # Define colors: we create a vector of color values, one for
  # each feature, and we give it names of the feature ID. Then we
  # can easily get the color value from the feature name.
  # A: make a vector of color values. The syntax may appear unusual -
  #    colorRampPalette() returns a function, and we simply append
  #    the parameter (number-of-features) without assigning the function
  #    to its own variable name.
  ftrCol <- colorRampPalette(c("#f2003c", "#F0A200", "#f0ea00",
                               "#62C923", "#0A9A9B", "#1958C3",
                               "#8000D3", "#D0007F"),
                             space="Lab",
                             interpolate="linear")(nrow(DB$feature))
  # B: Features may overlap, so we make the colors transparent by setting
  #    their "alpha channel" to 1/3  (hex: 55)
  ftrCol <- paste0(ftrCol, "55")
  # C: we asssign names
  names(ftrCol) <- DB$feature$ID
  # E.g. color for the third feature: ftrCol[ DB$feature$ID[3] ]

  # find the row-index of the protein ID in the protein table of DB
  iProtein <- which(DB$protein$name == name)

  # write the name of the protein
  text(-30, y, adj=1, labels=name, cex=0.75 )

  #draw a line from 0 to nchar(sequence-of-the-protein)
  lines(c(0, nchar(DB$protein$sequence[iProtein])), c(y, y),
        lwd=3, col="#999999")

  # get the rows of feature annotations for the protein
  iFtr <- which(DB$annotation$proteinID == DB$protein$ID[iProtein])

  # draw a colored box for each feature
  for (i in iFtr) {
    drawBox(DB$annotation$start[i],
            DB$annotation$end[i],
            y,
            ftrCol[ DB$annotation$featureID[i] ])
  }
}

# Plot each annotated protein:
# Get the rows of all unique annotated Mbp1 proteins in myDB

iRows <- grep("^MBP1_", myDB$protein$name)

# define the size of the plot-frame to accomodate all proteins
yMax <- length(iRows) * 1.1
xMax <- max(nchar(myDB$protein$sequence[iRows])) * 1.1  # longest sequence

# plot an empty frame
plot(1, 1,
     xlim = c(-200, xMax + 100),
     ylim = c(0, yMax),
     type = "n",
     axes = FALSE,
     bty = "n",
     main = "Mbp1 orthologue domain annotations",
     xlab = "sequence position",
     ylab="")
axis(1, at = seq(0, xMax, by = 100))
myCol <- colorRampPalette(c("#f2003c", "#F0A200",
                            "#f0ea00", "#62C923",
                            "#0A9A9B", "#1958C3",
                            "#8000D3", "#D0007F"),
                          space="Lab",
                          interpolate="linear")(nrow(myDB$feature))
myCol <- paste0(myCol, "55")
legend(xMax - 150, 6,
       legend = myDB$feature$name,
       cex = 0.7,
       fill = myCol)


# Finally, iterate over all proteins and call plotProtein()
for (i in seq_along(iRows)) {
  plotProtein(myDB, myDB$protein$name[iRows[i]], i)
}

# The plot shows what is variable and what is constant about the annotations in
# a group of related proteins. Your MBP1_MYSPE annotations should appear at the
# top.

# Task:
#    Put a copy of the plot into your journal and interpret it with respect
#    to MBP1_MYSPE, i.e. and note what you learn about MBP1_MYSPE from the plot.



# [END]
