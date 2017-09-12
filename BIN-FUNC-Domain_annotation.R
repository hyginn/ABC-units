# BIN-FUNC-Domain_annotation.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-FUNC-Domain_annotation unit.
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

# = 1 SMART Domain annotations

# Plot domain annotations as colored rectangles on a sequence.
# Step one: enter your domain annotations as features into the database.
#
# == Update myDB
# If the reference database has changed, we need to merge it in with myDB.
load("myDB.03.RData")  # load the previous version of myDB
# the new version of refDB was loaded when you
# pulled it from GitHub, and then typed init()
myDB <- dbMerge(myDB)  # merge the two databases and update myDB with the result
save(myDB, file = "myDB.04.RData") # save the new version

# == Update myDB

# Every annotated feature requires its own entry in the database. You have added
# the feature for the "APSES fold" before, so you can copy and edit that code
# from your myCode.R script. Here is again the table of feature IDs:
myDB$feature[ , c("ID", "name", "description")]

# Add every SMART annotated feaure for MBP1_YFO to the database. If you make
# mistakes, just reload the latest version (probably "myDB.04.RData"), then run
# your corrected annotation script again. Execute ...
myDB$proteinAnnotation
# ... to confirm.
#
# Once you are sure your annotations are correct, save the database again.
save(myDB, file = "myDB.05.RData") # save the new version
#
# Now let's plot the annotations.
#
# We need a small utility function that draws the annotation boxes on a
# representation of sequence. It will accept the left and right boundaries, the
# height and the color of the box and plot it using R's rect() function.

drawBox <- function(xLeft, xRight, y, colour) {
  # Draw a box from xLeft to xRight at y, filled with colour
  rect(xLeft, (y - 0.1), xRight, (y + 0.1),
       border = "black", col = colour)
}

# test this:
plot(c(-1.5, 1.5), c(0, 0), type = "l")
drawBox(-1, 1, 0.0, "peachpuff")

# Next, we define a function to plot annotations for one protein: the name of
# the protein, a horizontal grey line for its length, and all of its features.

plotProtein <- function(DB, ID, y) {
  # DB: protein database, probably you want myDB
  # ID: the ID of the protein to plot.
  # y: where to draw the plot
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
  #    their "alpha channel" to 1/2  (hex: 7F)
  ftrCol <- paste(ftrCol, "7F", sep = "")
  # C: we asssign names
  names(ftrCol) <- DB$feature$ID
  # E.g. color for the third feature: ftrCol[ DB$feature$ID[3] ]

  # find the row-index of the protein ID in the protein table of DB
  iProtein <- which(DB$protein$ID == ID)

  # write the name of the protein
  text(-30, y, adj=1, labels=DB$protein$name[iProtein], cex=0.75 )

  #draw a line from 0 to nchar(sequence-of-the-protein)
  lines(c(0, nchar(DB$protein$sequence[iProtein])), c(y, y),
        lwd=3, col="#999999")

  # get the rows of feature annotations for the protein
  iFtr <- which(DB$proteinAnnotation$protein.ID == ID)

  # draw a colored box for each feature
  for (i in iFtr) {
    drawBox(DB$proteinAnnotation$start[i],
            DB$proteinAnnotation$end[i],
            y,
            ftrCol[ DB$proteinAnnotation$feature.ID[i] ])
  }
}

# Plot each annotated protein:
# Get the rows of all unique annotated protein IDs in the protein table
iRows <- which(myDB$protein$ID %in% unique(myDB$proteinAnnotation$protein.ID))
# define the size of the plot-frame to accomodate all proteins
yMax <- length(iRows) * 1.1
xMax <- max(nchar(myDB$protein$sequence[iRows])) * 1.1  # longest sequence

# plot an empty frame
plot(1,1, xlim=c(-200, xMax), ylim=c(0, yMax),
     type="n", axes=FALSE, bty="n", xlab="sequence position", ylab="")
axis(1, at = seq(0, xMax, by = 100))

# Finally, iterate over all proteins and call plotProtein()
for (i in 1:length(iRows)) {
  plotProtein(myDB, myDB$protein$ID[iRows[i]], i)
}

# The plot shows clearly what is variable and what is constant about the
# annotations in a group of related proteins. Print the plot and bring it to
# class for the next quiz.
#

# = 1 Tasks




# [END]
