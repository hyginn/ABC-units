# tocID <- "BIN-FUNC-Domain_annotation.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-FUNC-Domain_annotation unit.
#
# ==============================================================================
# Version:  1.4
#
# Date:     2017-11  -  2020-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.4    Add code for shared data import from the Wiki
#           1.3    Add code for database export to JSON and instructions
#                  for uploading annotations to the Public Student Wiki page
#           1.2    Consistently: data in ./myScripts/ ;
#                    begin SHARING DATA section
#           1.1    2020 Updates
#           1.0    Live version 2017
#           0.1    First code copied from 2016 material.
#
# TODO:
#           Put the domain plot into a function
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
#TOC>   Section  Title                                                 Line
#TOC> ---------------------------------------------------------------------
#TOC>   1        Update your database script                             51
#TOC>   1.1        Preparing an annotation file ...                      58
#TOC>   1.1.1          BEFORE  "BIN-ALI-Optimal_sequence_alignment"      61
#TOC>   1.1.2          AFTER "BIN-ALI-Optimal_sequence_alignment"       109
#TOC>   1.2        Execute and Validate                                 136
#TOC>   2        Plot Annotations                                       161
#TOC>   3        SHARING DATA                                           287
#TOC>   3.1        Post MBP1_MYSPE as JSON data                         303
#TOC>   3.2        Import shared MBP1_MYSPE from the Wiki               326
#TOC> 
#TOC> ==========================================================================


# =    1  Update your database script  =========================================


# Since you have recorded domain features at the SMART database, we can store
# the feature annotations in myDB ...


# ==   1.1  Preparing an annotation file ...  ==================================


# ===   1.1.1  BEFORE  "BIN-ALI-Optimal_sequence_alignment"
#
#   IF YOU HAVE NOT YET COMPLETED THE BIN-ALI-OPTIMAL_SEQUENCE_ALIGNMENT UNIT:
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
#   - Save your file in the ./myScripts/ folder.
#
#   - Validate your file online at https://jsonlint.com/
#
#   - Update your "./myScripts/makeProteinDB.R" script to load your new
#     annotation when you recreate the database. Open the script in the
#     RStudio editor, and add the following command at the end:
#
#     myDB <- dbAddAnnotation(myDB,
#         jsonlite::fromJSON("./myScripts/<MYSPE>-Annotations.json"))
#                                         ^^^^^^^
#                                        edit this!
#
#   - save and close the file.
#
# Then SKIP the next section.
#
#
# ===   1.1.2  AFTER "BIN-ALI-Optimal_sequence_alignment"  
#
#   IF YOU HAVE ALREADY COMPLETED THE BIN-ALI-OPTIMAL_SEQUENCE_ALIGNMENT UNIT:
#
#   You SHOULD have a file called "<MYSPE>-Annotations.json" in the
#  ./myScripts/ directory:
#
#   - Open the file in the RStudio editor.
#
#   - Make as many copies of the "APSES fold" line as you have found
#     features in SMART.
#
#   - Add a comma after every line except for the last one
#
#   - Edit the annotations but include only features that are in the
#     myDB$feature table. Check which features are in the database by executing
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
#     source("./myScripts/makeProteinDB.R")
#
#     This should run without errors or warnings. If it doesn't work and you
#     can't figure out quickly what's happening, ask for help on the
#     Discussion Board.
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
# sequence, as an example of using the R plotting system for generic, data
# driven images.

# We need a small utility function that draws the annotation boxes on a
# representation of sequence. It should accept the start and end coordinates,
# the y value where it should be plotted and the color of the box, and plot a
# rectangle using R's rect() function.

drawBox <- function(xStart, xEnd, y, myCol, DELTA = 0.2) {
  # Draw a box from xStart to xEnd at y, filled with colour myCol
  # The height of the box is y +- DELTA
  rect(xStart, (y - DELTA), xEnd, (y + DELTA),
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
oPar <- par(mar = c(4.2, 0.1, 3, 0.1))  # save the current plot parameters and
                                        # decrease margins
plot(1, 1,
     xlim = c(-200, xMax + 100),
     ylim = c(0, yMax),
     type = "n",
     axes = FALSE,
     bty = "n",
     main = "Mbp1 orthologue domain annotations",
     xlab = "sequence position",
     cex.axis = 0.8,
     ylab="")
axis(1, at = seq(0, xMax, by = 100))
myCol <- colorRampPalette(c("#f2003c", "#F0A200",
                            "#f0ea00", "#62C923",
                            "#0A9A9B", "#1958C3",
                            "#8000D3", "#D0007F"),
                          space="Lab",
                          interpolate="linear")(nrow(myDB$feature))
myCol <- paste0(myCol, "55")
legend(xMax - 150, 7,
       legend = myDB$feature$name,
       cex = 0.7,
       fill = myCol,
       bty = "n")

# Finally, iterate over all proteins and call plotProtein()
for (i in seq_along(iRows)) {
  plotProtein(myDB, myDB$protein$name[iRows[i]], i)
}
par(oPar)  # reset the plot parameters


# The plot shows what is variable and what is constant about the annotations in
# a group of related proteins. Your MBP1_MYSPE annotations should appear at the
# top.

# Task:
#    Put a copy of the plot into your journal and interpret it with respect
#    to MBP1_MYSPE, i.e. and note what you learn about MBP1_MYSPE from the plot.

# Task:
#    It would be better to align the motif borders, at least approximately (not
#    all proteins have all motifs). How would you go about doing that?

# =    3  SHARING DATA  ========================================================

# It's particularly interesting to compare such annotations across many
# homologous proteins. I have created a page on the Student Wiki () that you can
# edit, and then download the data from the entire class directly to your
# RStudio project.
#

# I have provided a function that extracts all information that refers to a
# single protein from the database, and prints it out as well-formatted JSON,
# suitable to be pasted into our shareable Wiki-page. There is a fair amount of
# bookkeeping involved, but the code is not otherwise very enlightening so I
# will spare you the details - it's in "./scripts/ABC-dbUtilities.R" if you
# would want to have a look.


# ==   3.1  Post MBP1_MYSPE as JSON data  ======================================

# Task:
# =====
# 1: Run the following code:

cat("{{Vspace}}",
    "<!-- ==== BEGIN  PROTEIN ==== -->",
    "<pre class=\"protein-data\">",
    dbProt2JSON(sprintf("MBP1_%s", biCode(MYSPE))),
    "</pre>",
    "<!-- ===== END PROTEIN ====== -->",
    "", sep = "\n"
)

# 2: Copy the entire output from the console.
# 3: Navigate to
#      http://steipe.biochemistry.utoronto.ca/abc/students/index.php/Public
#    ... edit the page, and paste your output at the top.
# 4: Save your edits.



# ==   3.2  Import shared MBP1_MYSPE from the Wiki  ============================

# Once we have collected a number of protein annotations, we can access the
# Wiki-page and import the data into our database. The Wiki page is  an html
# document with lots of MediaWiki specific stuff - but the contents we are
# interested in is enclosed in <pre class="protein-data"> ... </pre> tags. These
# work like normal HTML <pre> tags, but we have defined a special class for them
# to make it easy to parse out the contents we want. The rvest:: package in
# combination with xml2:: provides us with all the tools we need for such
# "Webscraping" of data....

if (! requireNamespace("rvest", quietly=TRUE)) {
  install.packages("rvest")
}

if (! requireNamespace("xml2", quietly=TRUE)) {
  install.packages("xml2")
}

# Here's the process:
# The URL is an "open" page on the student Wiki. Users that are not logged in
# can view the contents, but you can only edit if you are logged in.
myURL <- "http://steipe.biochemistry.utoronto.ca/abc/students/index.php/Public"

# First thing is to retrieve the HTML from the url...
x <- xml2::read_html(myURL)

# This retrieves the page source, but that still needs to be parsed into its
# logical elements. HTML is a subset of XML and such documents are structured as
# trees, that have "nodes" which are demarcated with "tags". rvest::html_nodes()
# parses out the document structure and then uses a so-called "xpath" expression
# to select nodes we are interested in. Now, xpath is one of those specialized
# languages of which there are a few more to learn than one would care for. You
# MUST know how to format sprintf() expressions, and you SHOULD be competent
# with regular expressions. But if you want to be really competent in your work,
# basic HTML and CSS is required ... and enough knowledge about xpath to be able
# to search on Stackoverflow for what you need for parsing data out of Web
# documents...

# The expression we use below is:
#   - get any node anywhere in the tree ("//*") ...
#   - that has a particular attribute("[@ ... ]").
#   - The attribute we want is that the class of the node is "protein-data";
#      that is the class we have defined for our <pre> tags.
# As a result of this selection, we get a list of pointers to the document tree.
y <- rvest::html_nodes(x, xpath ='//*[@class="protein-data"]')

# Next we fetch the actual payload - the text - from the tree:
# rvest::html_text() gets the text from the list of pointers. The result is a
# normal list of character strings.
z <- rvest::html_text(y)

# Finally we can iterate over the list, and add all proteins we don't already
# have to our database. There may well be items that are rejected because they
# are already present in the database - for example, unless somebody has
# annotated new features, all of the features are already there. Don't worry -
# that is intended; we don't want duplicate entries.

for (thisJSON in z) {
  thisData <- jsonlite::fromJSON(thisJSON)
  if (! thisData$protein$name %in% myDB$protein$name) {
    myDB <- dbAddProtein(myDB, thisData$protein)
    myDB <- dbAddTaxonomy(myDB, thisData$taxonomy)
    myDB <- dbAddFeature(myDB, thisData$feature)
    myDB <- dbAddAnnotation(myDB, thisData$annotation)
  }
}

# Finally, we can repeat our domain plot with the results - which now includes the shared proteins:

iRows <- grep("^MBP1_", myDB$protein$name)
yMax <- length(iRows) * 1.1
xMax <- max(nchar(myDB$protein$sequence[iRows])) * 1.1  # longest sequence

# plot an empty frame
oPar <- par(mar = c(4.2, 0.1, 3, 0.1))
plot(1, 1,
     xlim = c(-200, xMax + 100),
     ylim = c(0, yMax),
     type = "n",
     axes = FALSE,
     bty = "n",
     main = "Mbp1 orthologue domain annotations",
     xlab = "sequence position",
     cex.axis = 0.8,
     ylab="")
axis(1, at = seq(0, xMax, by = 100))
myCol <- colorRampPalette(c("#f2003c", "#F0A200",
                            "#f0ea00", "#62C923",
                            "#0A9A9B", "#1958C3",
                            "#8000D3", "#D0007F"),
                          space="Lab",
                          interpolate="linear")(nrow(myDB$feature))
myCol <- paste0(myCol, "55")
legend(xMax - 150, 7,
       legend = myDB$feature$name,
       cex = 0.7,
       fill = myCol,
       bty = "n")

for (i in seq_along(iRows)) {
  plotProtein(myDB, myDB$protein$name[iRows[i]], i)
}
par(oPar)  # reset the plot parameters

# ... the more proteins we can compare, the more we learn about the
# architectural principles of this family's domains.


# [END]
