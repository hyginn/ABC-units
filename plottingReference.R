# tocID <- "plottingReference.R"
#
# Purpose: Reference to graphical output in R.
#
#
# Version: 2.0
#
# Date:    2017  09  -  2020 10
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 2.0    Comprehensive reference with basic and advanced options based on
#          an integrated yeast gene expression dataset. Full rewrite of most
#          sections.
# V 1.1    Stylistic improvements, code polish, and additional examples
# V 1.0    Digest from "plotting reference" files
#
# ToDo:
#
#


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                              Line
#TOC> ------------------------------------------------------------------
#TOC>   01       INITIALIZE                                           67
#TOC>   02       THIS REFERENCE ...                                   72
#TOC>   02.1       Dataset Documentation                              79
#TOC>   03       PROPORTIONS AND DISTRIBUTIONS                       191
#TOC>   03.1       barplot()                                         196
#TOC>   03.2       pie()                                             220
#TOC>   03.3       boxplot()                                         252
#TOC>   03.4       hist()                                            310
#TOC>   03.4.1         overlaying histograms                         355
#TOC>   04       THE plot() FUNCTION                                 398
#TOC>   04.1       line plots                                        402
#TOC>   05       ENCODING INFORMATION: SYMBOL, SIZE, COLOR           501
#TOC>   05.1       pch ("plotting character" symbols)                507
#TOC>   05.1.1         Line types                                    682
#TOC>   06       COLOUR                                              704
#TOC>   06.1       Colours by number                                 712
#TOC>   06.2       Colours by name                                   730
#TOC>   06.3       Colours as hex-triplets                           756
#TOC>   06.3.1         Inbuilt palettes                              811
#TOC>   06.3.2         Custom palettes                               887
#TOC>   06.3.3         Transparency: The Alpha Channel               937
#TOC>   06.4       abline(), lines()  and segments()                 983
#TOC>   07       AXES                                               1017
#TOC>   08       LEGENDS                                            1054
#TOC>   08.1       basic legends                                    1057
#TOC>   08.2       Color bars                                       1061
#TOC>   09       LAYOUT                                             1180
#TOC>   10       TEXT                                               1215
#TOC>   11       DRAWING ON PLOTS                                   1243
#TOC>   12       IMAGES                                             1310
#TOC>   13       CONTOUR LINES                                      1315
#TOC>   14       3D PLOTS                                           1320
#TOC>   15       GRAPHS AND NETWORKS                                1325
#TOC>   16       OTHER GRPAHICS PACKAGES                            1330
#TOC>   17       INTERACTIVE PLOTS                                  1354
#TOC>   17.1       locator()                                        1358
#TOC>   17.2       plotly::                                         1361
#TOC>
#TOC> ==========================================================================


# =    01  INITIALIZE  =========================================================

SC <- readRDS("data/SC.rds")   # <<<- execute this line first


# =    02  THIS REFERENCE ...  =================================================

# This script covers basic plotting and graphical output options of R. The
# functions are demonstrated with an integrated data set of yeast gene
# expression data and annotations.


# ==   02.1  Dataset Documentation  ============================================

# You do not need to study the dataset in detail but it will be useful to refer
# to this documentation from time to time to understand what data is being
# plotted and what one can therefore learn from the results.

# The integrated dataset SC is a list that contains an analysis of yeast gene
# expression profiles originally published by Pramila et al. (2002; PMID:
# 12464633) and accessible as GSE3653 on GEO
# (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3635). This is a
# high-resolution (5 min. interval) data set spanning two cycles of a
# synchronized yeast culture. It includes 6228 expression profiles with 25
# time-points assigned to yeast systematic names. It is complemented with gene
# annotations taken from SGD, and augmented with a set of GO slims. GO slims
# were rerieved by navigating to https://www.yeastgenome.org/ and using the menu
# to choose Function >> Gene Ontology >> GOslim Mapping File. This downloads
# go_slim_mapping.tab which was further processed to annotate each gene with its
# corresponding GO terms. To evaluate cyclical regulation, a curve-fitting
# protocol was followed to fit a cyclical epxpression model with parameters
# amplitude, phase, frequency, exponential damping, and baseline shift. After
# fitting the profiles, correlations with the fitted model were calculated and
# the first peak of the model was determined as a marker of when the expression
# was ON in the cell-cycle. Finally, the expression peak markers were assigned
# to seven phases of the cell cycle, and GO term enrichments were computed. A
# A subset was selected with the following parameters:

# selCCl <-
#   nlsParams$A > 0.075 &   # reasonably high Amplitude
#   nlsParams$cor > 0.6 &   # good correlation with the parametrized model
#   nlsParams$f > 0.75 &    # period between 0.75 ...
#   nlsParams$f < 1.333 &   # ... and 1.333 hours
#   nlsParams$k < 0.03 &    # limiting the exponential damping
#   nlsParams$k > -0.001
# sum(selCCl)               # 1,297 genes in the data set show some level
#                           # of periodic expression-variation in the cell
#                           # cycle.
#
# Thus the dataset comprises numeric and categorical data
# Here are the data details:
# SC$xpr     :  expression profiles. Numeric matrix with 1,297 rows of genes
# ======                             and 25 colums of time points in 5 minute
#                                    intervals
# Example:
# --------
# SC$xpr[169, 1:5]
# >          t.0          t.5         t.10         t.15         t.20
# > -0.055703429 -0.116124955 -0.135965733 -0.068665035 -0.001819363
#
#
# SC$mdl     :  parametrized model. Data frame with 1,297 rows of genes
# ======
# SC$mdl$A     : model parameter: amplitude   (log ratio)
# SC$mdl$phi   : model parameter: phase       (degrees)
# SC$mdl$f     : model parameter: 1/frequency (hours)
# SC$mdl$k     : model parameter: damping
# SC$mdl$B     : model parameter: baseline shift
# SC$mdl$cor   : correlation of observation and model
# SC$mdl$peaks : timepoint of first expression peak of model (minutes)
#
# # Example:
# --------
# SC$mdl[169, ]
# >         A       phi        f          k           B       cor    peaks
# > 0.1007908 -81.04413 1.012092 0.00656052 -0.01857916 0.8437175 37.52217
#
#
# SC$ann     :  annotations.  Data frame with 1,297 rows of genes
# ======
# SC$ann$SGD          :  Saccharomyces Gene Database identifier
# SC$ann$sysName      :  Yeast gene systematic name
# SC$ann$stdName      :  The standard name under which the gene is known
# SC$ann$alias        :  common alias name(s)
# SC$ann$description  :  free-text description of the gene
# SC$ann$GO           :  GO IDs of the SGD GO slim subset annotated to the gene
#
# Example:
# --------
# SC$ann[169, ]
# >        SGD sysName stdName                     alias
# > S000002214 YDL056W    MBP1 transcription factor MBP1
# > description
# > YDL056W Transcription factor; involved in regulation [...]
# > GO
# > YDL056W GO:0005634 GO:0003677 GO:0001071 GO:0000278
#
# SC$phases
# ---------
# Data frame with a definition of yeast cell cycle phase names and
# time-points derived from inspection of the  SC$mdl$peaks  distribution.
#         names start ends
#   1     Sense     0    6
#   2      Prep     6   15
#   3 Replicate    15   30
#   4  Assemble    30   37
#   5 Segregate    37   43
#   6 Duplicate    43   55
#   7 Stabilize    55   65
#
# SC$GO
# -----
# Data frame with GO term information and GO term enrichment data
# SC$GO$ontology  : {C|F|P} identifying the "cellular component", "molecular
#                           function" or "biological process" ontology
# SC$GO$GOid      : Gene Ontyology ID of the term
# SC$GO$label     : Short label of the term
# SC$GO$tAll      : Count of times the term is annotated in any gene
# SC$GO$tCC       : ... times term is annotated to a gene in this dataset
# SC$GO$t...      : ... annotated to one ofg the cell-cycle phases
# SC$GO$xs...     : ... enrichment factor in each cell cycle phase
#


# =    03  PROPORTIONS AND DISTRIBUTIONS  ======================================

# Distributions of numeric variable characterize the values. Proportions compare values. A typical use case is to characterize a set of measurements, or compare several sets of measurements.


# ==   03.1  barplot()  ========================================================
# Draws a bar with height proportional to a given number.

barplot(mean(SC$xpr[ , "t.20"])) # barplot of the mean expression at t.20

# Barplots show only a single number. They tell us nothing about the underlying distribution. They are commonly used to compare several distributions:

barplot(colMeans(SC$xpr[ , c("t.0", "t.20", "t.40", "t.60")]))

# barplots take the usual plotting parameters
myPal <- colorRampPalette(c("#00FF0055", "#0066FF55",    # a "palette" function
                            "#00FFFF55", "#CCCCAA55",
                            "#FF00FF55", "#FF006655",
                            "#FFFFFF55"), alpha = TRUE)  # transparent colors

barplot(colMeans(SC$xpr[ , 1:13]),         # mean expression changes, first hour
        main = "Mean expression changes",
        xlab = "time points",
        cex.names = 0.6,                   # scale column names
        ylab = "Expression (log ratio)",
        cex.axis = 0.8,                    # scale of y-axis labels
        col = myPal(13))                   # get 13 color values from myPal()


# ==   03.2  pie()  ============================================================
# Generally not preferred, but simple to do

table(SC$GO$ontology)
# Example: how many genes are anotated to the GO terms in the Cellular Component
# category?
CCids <- SC$GO$GOid[SC$GO$ontology == "C"]  # get GO ids for "C" ontology
x <- unlist(strsplit(SC$ann$GO, " ")) # split GO term strings on blank-spaces
                                      # and dump them all into one single vector
x <- x[x %in% CCids] # subset the ones in CCIDs (about 40%)

pie(table(x), cex = 0.5)

# This can be better visualized in a barplot
myT <- sort(table(x), decreasing = TRUE)
b <- barplot(myT,                  # assign the plot - we need the x-coordinates
     main = "Genes annotated to GO terms in the Cellular Component Ontology",
     cex.main = 0.8,
     names.arg = "",               # blank the names - otherwise they are
     ylab = "counts",              #   taken from the names of the object
     cex.axis = 0.8,
     col = colorRampPalette(c("#CCAAFF", "#FFDDFF", "#FFFFFF"))(length(myT)),
     ylim = c(0, 1.5 * max(myT)))  # make  space at the top

# Add the GOids as text()
text(b, myT,                 # b holds the x coordinates, myT is y
     labels = names(myT),    # take the names from the table - names
     srt = 90,               # rotate to vertical
     adj = c(-0.1, 0.5),     # align center, top
     cex = 0.5)


# ==   03.3  boxplot()  ========================================================
# A boxplot() is almost always preferred to a barplot, since it includes an
# estimate of the distribution:

boxplot(SC$xpr[ , "t.20"])   # boxplot of expression values at t.20

# What are these elements?
# Assigning the plot to a variable makes the numbers available so we can
# use them for annotation:
oPar <- par(mar = c(0.2,3,0.2,0.2))   # reduce the margins
( b <- boxplot(SC$xpr[ , "t.20"]))

xT <- 1.3
cT <- "#DD0000"
# outliers
myOut    <- mean(b$out[b$out > b$stats[5]]) # outliers above the whiskers
myOut[2] <- mean(b$out[b$out < b$stats[1]]) # outliers below the whiskers
text(1.04, myOut,        adj=c(0,0.5), cex=0.7, col=cT, "outliers")

# whiskers
text(1.12, b$stats[1,1], adj=c(0,0.5), cex=0.7, col=cT, "min(x_i > Q1-1.5*IQR)")
text(1.12, b$stats[5,1], adj=c(0,0.5), cex=0.7, col=cT, "max(x_i < Q3+1.5*IQR)")

# quartiles
# "min(x_i > Q1-1.5*IQR)" is "the smallest x in the data that still lies
# within 1.5 times the inter-quartile range below Q1
text(1.23, b$stats[2,1], adj=c(0,0.5), cex=0.7, col=cT, "Q1 (lower quartile)")
text(1.23, b$stats[4,1], adj=c(0,0.5), cex=0.7, col=cT, "Q3 (upper quartile)")

# median
text(1.23, b$stats[3,1], adj=c(0,0.5), cex=0.9, col=cT, "median")

# the mean is not shown
text(0.79, mean(SC$xpr[,"t.20"]), adj=c(1,0.5), cex=0.7, col="#999999", "mean")

par(oPar)  # reset the margins

# If we have datasets with more columns, each column gets its own boxplot
boxplot(SC$xpr[ , c("t.0", "t.20", "t.40", "t.60")])

# boxplots take the usual plotting parameters
myPal <- colorRampPalette(c("#00FF0055", "#0066FF55",    # a "palette" function
                            "#00FFFF55", "#CCCCAA55",
                            "#FF00FF55", "#FF006655",
                            "#FFFFFF55"), alpha = TRUE)  # transparent colors

boxplot(SC$xpr[ , 1:13],                   #  expression changes, first hour
        main = "Expression levels",
        xlab = "time points",
        cex.axis = 0.6,                    # scale axis names
        ylab = "Expression (log ratio)",
        col = myPal(13))                   # get 13 color values from myPal()

# This plot shows that there is a noticable fluctuation of expression levels
# over the course of the first 60 minutes of the experiment. This may be due to
# these genes being selected as having cyclically varying expression profiles.


# ==   03.4  hist()  ===========================================================
# Histograms are superbly informative and one of the workhorses of our analyses. A histogram answers the question how many observations do we have in a particular range of a (continuously varying) variable.

# Here we plot a histogram of the first expression peaks of genes, computed as
# the first peak of the fitted model for every gene in our set. We might be
# interested in whether the expression of different genes peaks continuously
# over the cycle, or whetehr there are time intervals where waves of expression
# of different genes can be observed in response to cycle-control signals.
hist(SC$mdl$peaks)

# We can explicitly set breakpoints for the histogram:
# Here we set breaks at 4 minute intervals
( myBreaks <- seq(0, 80, by = 4) )  # Note: 21 values bound 20 intervals
hist(SC$mdl$peaks, breaks = myBreaks)

# Assigning the output of hist() makes the values
# used in constructing the histogram accessible:

( H <- hist(SC$mdl$peaks) )

# for example, we can use this to plot the actual numbers:
text(H$mids,                       # midpoints of the bars
     H$counts,                     # counts -> probabilities
     labels = H$counts,            # add labels to each bar: count numbers
     adj = c(0.5, -0.5),           # text is centered on x and raised in y
     cex = 0.6,
     col = "#8888FF")

# histograms take the usual plotting parameters
myPal <- colorRampPalette(c("#00FF0055", "#0066FF55",    # a "palette" function
                            "#00FFFF55", "#CCCCAA55",
                            "#FF00FF55", "#FF006655",
                            "#FFFFFF55"), alpha = TRUE)  # transparent colors

myBreaks <- seq(0, 80, length.out = 81)    # one-minute intervals
hist(SC$mdl$peaks,                         # gene expression first peaks
     breaks = myBreaks,
     main = "Timing of first peaks of during the cell-cycle",
     sub = "The bimodal distribution distinguishes replication from division",
     cex.sub = 0.8,
     xlab = "t (min.)",
     ylab = "counts",
     col = myPal(81))                      # get 81 color values from myPal()


# ===   03.4.1  overlaying histograms

# Histograms can be plotted one-over another to compare them
# Example: when do genes with different GO term annotations first peak?

# GO:0005730 "nucleolus"
# GO:0032200 "telomere organization"
# GO:0006260 "DNA replication"
sel1 <- grep("GO:0005730", SC$ann$GO)  # 103 genes
sel2 <- grep("GO:0032200", SC$ann$GO)  # 32 genes
sel3 <- grep("GO:0006260", SC$ann$GO)  # 67 genes

# define breaks. This is important to make the histograms comparable
myBreaks <- seq(0, 80, by = 4)
# plot the first histogram

hist(SC$mdl$peaks[sel1], breaks = myBreaks,
     main = "Peak timing for different GO terms",
     xlab = "t (min)",
     ylab = "Counts",
     ylim = c(0, 50),
     col = "#FF007955")

# first overlay
hist(SC$mdl$peaks[sel2], breaks = myBreaks,
     col = "#00B2FF55",
     add = TRUE)

# second overlay
hist(SC$mdl$peaks[sel3], breaks = myBreaks,
     col = "#7FDFC955",
     add = TRUE)

# Legend
legend("topright",
       fill = c("#FF007955", "#00B2FF55", "#7FDFC955"),
       legend = c("GO:0005730\n(nucleolus)\n",   # "\n" is a line-break
                  "GO:0032200\n(telomere organization)\n",
                  "GO:0006260\n(DNA replication)\n"),
       cex = 0.6,
       bty = "n")  # no box around the legend


# =    04  THE plot() FUNCTION  ================================================
# plot() is the workhorse of data visualization in R. But plot() is also a "generic" - i.e. different "classes" can define their own plot-methods and plot will "dispatch" the right method for a class.


# ==   04.1  line plots  =======================================================

# Lines are useful to visually associate ralated data points, like points in a
# sequence of events.

# plot an expression profile: value over time
plot(SC$xpr[169, ])

# note the x is implied and plotted as "index", the expression values are used
# for the y-axis by default.

# plot types
plot(SC$xpr[169, ], type = "p")  # points, the default
plot(SC$xpr[169, ], type = "l")  # lines
plot(SC$xpr[169, ], type = "b")  # both points and lines
plot(SC$xpr[169, ], type = "c")  # empty points and lines
plot(SC$xpr[169, ], type = "o")  # overplotted
plot(SC$xpr[169, ], type = "s")  # steps
plot(SC$xpr[169, ], type = "h")  # like histogram
plot(SC$xpr[169, ], type = "n")  # none - useful to create an empty frame
                                 # to add other elements later
text(SC$xpr[169, ], colnames(SC$xpr), cex = 0.5) # ... such as labels

# title and axis labels
iRow <- 169
t <- seq(0, 120, by = 5)               # set actual x-values
plot(t, SC$xpr[iRow, ],
     type = "b",
     main = "Expression profile",
     sub = "Data: GSE3653 (Pramila et al. 2002)",
     xlab = "time(min)",
     ylab = "Expression (log ratio)")

# adjusting font sizes
plot(t, SC$xpr[iRow, ],
     type = "b",
     main = "Expression profile",
     sub = "Data: GSE3653 (Pramila et al. 2002)",
     xlab = "time(min)",
     ylab = "Expression (log ratio)",
     cex.main = 1.1,                                # title
     cex.sub = 0.6,                                 # subtitle
     cex.axis = 0.7,                                # axis values
     cex.lab = 0.8)                                 # axis labels

# x and y axis ranges can be defined with a two-element vector
plot(t, SC$xpr[iRow, ],
     type = "b",
     xlim = c(0, 120),                      # tight against the values
     ylim = c(-max(abs(SC$xpr[iRow, ])),
               max(abs(SC$xpr[iRow, ]))),  # symmetric around 0
     main = "", sub = "", xlab = "", ylab = "")
abline(h = 0, col = "#FF000044")

# pull in annotations with sprintf()
plot(t, SC$xpr[iRow, ],
     type = "b",
     main = sprintf("Expression profile: %s (%s)",
                    SC$ann$stdName[iRow],
                    SC$ann$sysName[iRow]),
     sub = "", xlab = "", ylab = "")

# When exploring data, it is useful to have a standard view of data as a
# function:

xprPlot <- function(iRow) {
  t <- seq(0, 120, by=5)
  yMax <- max(abs(SC$xpr[iRow, ])) * 1.1

  plot(t, SC$xpr[iRow, ],
       ylim = c(-yMax, yMax),
       type = "b",
       main = sprintf("Expression profile: %s (%s)",
                      SC$ann$stdName[iRow],
                      SC$ann$sysName[iRow]),
       xlab = "time(min)",
       ylab = "Expression (log ratio)",
       cex.main = 1.1,
       cex.axis = 0.7,
       cex.lab = 0.8)
  abline(h = 0, col = "#FF000044")
}

xprPlot(169)
xprPlot(1124)

# better to visualize by plotting one OVER the other.
# This is done with the function points()

xprPlot(1124)
points(t, SC$xpr[169, ], type = "b", cex = 0.7, col = "#EE5500CC")



# In general, when we plot(), we produce a "scatterplot" of X and y values.
# scatterplot: one against the other
plot(SC$xpr[169, ], SC$xpr[1124, ])  # these two profiles are anticorrelated


# =    05  ENCODING INFORMATION: SYMBOL, SIZE, COLOuR  =========================

# We have many options to encode information in scatter plots (and others), in order to be able to visualize and discover patterns and relationships. The most important ones are the type of symbol to use, its size, and its color.



# ==   05.1  pch ("plotting character" symbols)  ===============================
# There are many types of symbols available to plot points and they can all be varied by type, size, and color to present additional information. Which one to use is defined in the argument  pch  . The default is pch = 1: the open circle. Characters 1:20 are regular symbols:

# reduce margins
oPar <- par(mar = c(1, 1, 1, 1))

# Empty plot frame ...
plot(c(0,11), c(0,11), type = "n", axes = FALSE, xlab = "", ylab = "")

# coordinates for first 25 symbols
x1 <- c(1:10, 1:10)
y1 <- c(rep(10, 10),rep(9, 10))
myPch <- 1:20

points(x1, y1, pch = myPch, cex = 1.2)

# pch 21:25 can have different border and fill colours. Borders are defined with
# "col", fills are defined with "fill":

x2 <- 1:5
y2 <- rep(8,5)

myCols <- c("#4B0055", "#00588B", "#009B95", "#53CC67", "#FDE333")

points(x2, y2, pch=21:25, col="#708090", bg = myCols, cex = 1.8)

# ten more symbols are defined as characters
x3 <- 1:10
y3 <- rep(6,10)
myPch <- c(".", "o", "O", "0","a","A", "*", "+","-","|")
points(x3, y3, pch = myPch) # note: "myPCh" is a character vector while
                            # the other pch were identified by integers

# The ASCII codes for characters 32 to 126 can also be used as plotting symbols
myPch <- 32:126
x4 <- (((myPch - 32) %% 19) + 2) / 2
y4 <- 5 - ((myPch - 32) %/% 19)
points(x4, y4, pch=32:126, col="#0055AA")

# I don't think I have ever actually used characters as plotting characters in
# this way; it is much more straightforward to just use the text() function in
# that case, rather than having to look up which ASCII code is which character,
# and remembering which ones are available as numbers, and which ones are
# identified by the character itself. But the first twenty five I use a lot, and
# I also end up looking up their codes a lot. So here is a function to plot them for reference. You can paste it into your .Rprofile to keep it at hand:

pchRef <- function() {
  # a reference plot for the first twenty five plotting characters

  oPar <- par(mar = c(1, 1, 1, 1), mfrow = c(1, 1))
  plot(c(0,6), c(0,6), type = "n", axes = FALSE, xlab = "", ylab = "")

  x <- ((0:24) %% 5) + 1
  y <- rep(5:1, each = 5)

  idx <- 1:20
  points(x[idx], y[idx], pch = idx, cex = 2)

  idx <- 21:25
  points(x[idx], y[idx], pch = idx, cex = 2, col = "#708090",
         bg=c("#CC0033AA", "#991566AA", "#662A99AA", "#323FCCAA", "#0055FFAA"))

  text(x - 0.3, y, labels = sprintf("%d:", 1:25), col = "#7766AA", cex = 0.8)

  par(oPar)
}


# reset margins
par(oPar)


# In addition to the default fonts, a large set of Hershey vector fonts is
# available which gives access to many more plotting and labeling options via
# text()
demo(Hershey)

# Plotting other symbols:
# In the most general way, Unicode characters can be plotted as text.
# The code is passed in hexadecimal, long integer, with a negative sign.
# Here is a quarter note (Unicode: 266a) using plot()
plot(0.5,0.5, pch=-0x266aL, cex=5, xlab="", ylab="")

# However, rendering may vary across platforms since it depends on unicode
# support. If your character can't be rendered, you'll only get an empty box.


# 1.1.1 Using pch to emphasize different categories

# In general, when we plot(), we produce a "scatterplot" of X and y values.
# scatterplot: one against the other

tRep <- c("t.15", "t.20", "t.25", "t.30")   # Replication-phase time points
tDup <- c("t.40", "t.45", "t.50", "t.55")   # Duplication-phase time points

plot(rowMeans(SC$xpr[ , tRep]), rowMeans(SC$xpr[ , tDup]))

# select top ten replication enriched GO terms
sel <- order(SC$GO$xsReplicate, decreasing = TRUE)[1:10]
SC$GO$label[sel]  # what are these?
GOrep <- SC$GO$GOid[sel]   # GO terms that are enriched in the Replication phase

sel <- order(SC$GO$xsDuplicate, decreasing = TRUE)[1:10]
SC$GO$label[sel]  # what are these?
GOdup <- SC$GO$GOid[sel]   # GO terms that are enriched in the Duplication phase

# define sets of genes that are annotated to one or the other GO id set
repGenes <- logical(nrow(SC$ann))             # prepare an empty vector
for (id in GOrep) {                           # whenever we find one of the
  repGenes <- repGenes | grepl(id, SC$ann$GO) # characteristic GO IDs annotated
}                                             # to a a gene, we put TRUE in its
sum(repGenes)                                 # spot of the repGenes vector.
xRep <- rowMeans(SC$xpr[ repGenes, tRep]) # repGenes in the replication phase
yRep <- rowMeans(SC$xpr[ repGenes, tDup]) # repGenes in the duplication phase

# same thing for the genes that are enriched in the duplication phase
dupGenes <- logical(nrow(SC$ann))
for (id in GOdup) {
  dupGenes <- dupGenes | grepl(id, SC$ann$GO)
}
sum(dupGenes)
xDup <- rowMeans(SC$xpr[ dupGenes, tRep]) # dupGenes in the replication phase
yDup <- rowMeans(SC$xpr[ dupGenes, tDup]) # dupGenes in the duplication phase

# Ready to plot:
# base plot without the ref... and dup... genes
basePlot <- function(...) {
  sel <- ! (repGenes | dupGenes)
  m <- sprintf("Expression changes between phases (%d genes)",
               nrow(SC$xpr))
  plot(rowMeans(SC$xpr[ sel, tRep]),
       rowMeans(SC$xpr[ sel, tDup]),
       main = m,
       xlab = "Expression in Replication phase",
       ylab = "Expression in Duplication phase",
       cex.main = 0.9,
       cex.lab = 0.8,
       ...)
  abline(h = 0, col = "#0088FF44")       # horizontal zero
  abline(v = 0, col = "#0088FF44")       # vertical zero
}

basePlot()
# plotting the rep.. and dup... sets with different pch
points(xRep, yRep, pch=24, col="#000000", bg="#FFFFFF")  # white up triangle
points(xDup, yDup, pch=25, col="#000000", bg="#FFFFFF")  # white down triangle

# These are clearly different distributions,  but the plot is overall too busy.

# Using color improves the plot:
basePlot(pch=21, col="#AAAAAA", bg="#DDDDDD")
points(xRep, yRep, pch=24, col="#000000", bg="#66DDFF")
points(xDup, yDup, pch=25, col="#000000", bg="#FF77FF")

# Even better may be to use transparent colors, to address the overlap:
basePlot(pch=19, col="#00000022")
points(xRep, yRep, pch=19, col="#33BBFF66")
points(xDup, yDup, pch=19, col="#FF44AA44")

# this needs a legend:
legend ("bottom",
        legend = c(
  sprintf("%d Genes with Replication-\nenriched GO terms\n ", sum(repGenes)),
  sprintf("%d Genes with Duplication-\nenriched GO terms\n ", sum(dupGenes)),
          "Other Genes\n "
          ),
        horiz = TRUE,
        pch = 19,
        col = c("#33BBFF66", "#FF44AA44", "#00000022"),
        cex = 0.6,
        pt.cex = 1.0,
        inset = 0.02,
        box.col = "#DDDDDD",
        bg = "#FFFFFF")

# ===   05.1.1  Line types

# Basically all plots take arguments lty to define the line type, and lwd
# to define line width

# empty plot ...
plot(c(0,10), c(0,10), type = "n", axes = FALSE, xlab = "", ylab = "")

# Line type
for (i in 1:8) {
  y <- 10.5-(i/2)
  segments(1,y,5,y, lty=i)
  text(6, y, paste("lty = ", i), col="grey60", adj=0, cex=0.75)
}

# Line width
for (i in 1:10) {
  y <- 5.5-(i/2)
  segments(1,y,5,y, lwd=(0.3*i)^2)
  text(6, y, paste("lwd = ", (0.3*i)^2), col="grey60", adj=0, cex=0.75)
}



# ==   05.1  cex ("character expansion" size)  =============================

# The size of characters can be controlled with the cex parameter. cex takes a
# vector of numbers, mapped to the vecors of plotted elements. The usual R
# conventions of vector recycling apply, so if cex s only a single number, it is
# applied to all plotted elements. Here is an example plotting expression
# amplitudes in a comparison of pre-replication vs. replication phase
# expression.


tPre <- c("t.0",  "t.5",  "t.10")   # Pre replication-phase time points
tMid <- c("t.25", "t.30", "t.35")   # Midway between replication and duplication

x <- rowMeans(SC$xpr[ , tPre])      # row means average out fluctuations
y <- rowMeans(SC$xpr[ , tMid])
plot(x, y)

myAmps <- SC$mdl$A                  # fetch the correlations of the fitted
                                    # low correlations are not well modeled
                                    # with a periodic function

hist(myAmps, breaks = 40)           # examine the distribution
hist(log(myAmps), breaks = 40)      # sometimes log() is be better suited
                                    # to map to symbol size

# useful cex sizes are approximately between 0.5 and 5
# Example:
N <- 10
cexMin <- 0.5
cexMax <- 5.0
myCex <- seq(cexMax, cexMin, length.out=N)  # Vector of cex values of length N
# We plot these points from largest to smallest, so we don't overplot
x1 <- runif(N)
y1 <- runif(N)
plot(x1, y1, xlim = 0:1, ylim = 0:1, xlab = "", ylab = "", pch = 16,
     cex = myCex, col= colorRampPalette(c("#F4F0F8", "#3355DD"))(N))

# Note that the smallest cex at 0.5 is _really_ small.
points(x1[N], y1[N], cex = 3.0, col = "#FF000077") # Here it is!

# To rescale a vector into the desired interval (lo, up) we write a little
# function. First we normalize the vector into the interval (0, 1), then we
# multiply it by up-lo, and add lo.
#

reScale <- function(x, lo = 0, up = 1) {
  return((((x-min(x)) / (max(x)-min(x))) * (up-lo)) + lo  )
}

# with this, our coding amplitudes as plot character size becomes quite easy:
plot(x, y, pch=16,
     cex = reScale(myAmps, lo=0.5, up=5),  # <<<---- That's all we need
     col="#4499BB22")
abline(a = 0, b = 1, col = "#DD333344")
abline(h = 0, col = "#33333344")
abline(v = 0, col = "#33333344")

# ToDo: add titles, scale legend, and interpretation


# =    06  COLOUR  =============================================================

# Colour is immensely important for data visualization,  and colours in R can be
# specified by number, by name, as hex-triplets, as rgb or hsv values, and
# through several color-space specifications. The can be used individually, as
# vectors that are defined manually, and through colour palettes. And the can be
# solid and transparent.

# ==   06.1  Colours by number  ================================================
# The col=... parameter for plots is 1 by default and you can
# set it to the range 0:8.
# 0: white
# 1: black (the default)
# 2: red
# 3: green
# 4: blue
# 5: cyan
# 6: magenta
# 7: yellow
# 8: grey
barplot(rep(1,9), col=0:8, axes=FALSE, names.arg=c(0:8))

# As you can see, the default primary colours are garish and offend even the
# most rudimentary sense of aesthetics. Using these colors won't earn you
# respect. Fortunately we have many more sophisticated ways to define colours.

# ==   06.2  Colours by name  ==================================================
# You may have noticed that "red", "green", and "blue" work for the col=...
# parameter, but you probably would not have imagined that "peachpuff",
# "firebrick" and "goldenrod" are valid as well. In fact, there are 657 named
# colours in R. Access them all by typing:
colors()

myCols <- c(
  "firebrick2",
  "tomato",
  "goldenrod1",
  "peachpuff",
  "papayawhip",
  "seashell",
  "whitesmoke"
)
pie(c(1, 1, 2, 3, 5, 8, 13),
    col = myCols,
    labels = myCols)

# Read more about named colours (and related topics) in Melanie Frazier's pdf:
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf

# I almost never use named colours anymore but use hex-triplets instead...


# ==   06.3  Colours as hex-triplets  ==========================================
# Hex triplets in R work exactly as in HTML: a triplet of
# RGB values in two-digit hexadecimal representation. The
# first two digits specify the red value, the second two
# are for green, then blue. R accepts a fourth pair of
# digits to optionally specify the transparency - the "Alpha" channel. The
# semantics of the hex-triplet is thus "#RRGGBB" or "#RRGGBBAA".
# Read more e.g. at http://en.wikipedia.org/wiki/Web_colors

# The function col2rgb() converts colour names to rgb values ...
col2rgb("violetred")

# ... and rgb() converts rgb values to hex-code:
rgb(1, 0.5, 0.23)

# to get hexcodes from names takes a bit more: divide and transpose -
rgb(t(col2rgb("violetred")/255))

# Here is a utility function to convert all manners of color-vectors to
# hex-triplets. You can copy it into your .Rprofile if you want  to keep it at hand.

col2hex <- function(colV, dput = FALSE, pal = FALSE, N = NULL) {
  # colV   a vector or scalar of color names, hex-codes or integers
  # dput   logical. If TRUE, character output is suitable for pasting into code.
  # pal    logical. If TRUE, the colors will be used in a palette. N must
  #                   be defined in that case
  # value: a vector of hex-codes of the same length, or a dput() string

  h <- rgb(t(col2rgb(colV) / 255))
  if (pal == TRUE) {
    if (is.null(N)) { stop("N must be defined if pal is TRUE.") }
    h <- colorRampPalette(h)(N)
  }
  if (dput == TRUE) {
    dput(h)
    return(invisible(h))
  }
  return(h)
}

# Examples:
myCols <- c( "firebrick2", "tomato",   "goldenrod1", "peachpuff",
             "papayawhip", "seashell", "whitesmoke")

col2hex(myCols)
col2hex(myCols, dput = TRUE)

# barplot with the original colors
barplot(rep(1, length(myCols)), col=myCols)

# expanding the palette with col2hex()
N <- 40
barplot(rep(1, N), col = col2hex(myCols,pal=TRUE,N=N) )


# ===   06.3.1  Inbuilt palettes

# In R, a palette is a function (!) that takes a number as its argument and
# returns that number of colors, those colors can then be used to color points
# in a plot or other elements. In principle, such colors can serve three
# distinct puposes:

# Qualitative colors have a different color for different types of entities.
# Here the color serves simply to distinguish the types, it represents
# "categorical information". For example, we can assign different colors to each
# individual amino acid.

# Sequential colors distinguish order in a sequence: the colors distinguish the first from the last element in the order. For example we often color the N-terminus of a protein blue(ish) and the C-terminus red(dish).

# Diverging colors encode some property and are scaled to represent particular
# values of the property. For example we might represent temperature with a
# black-body radiation color palette. Note that "sequential" can be seen as
# "diverging on order".

# R has are a number of inbuilt palettes; Here is a convenience function to view
# palettes. Input can be either a palette, the name of a hcl.colors palette, or a vector of colors:
plotPal <- function(pal, N = 40) {
  if (is.character(pal)) {
    if (length(pal) == 1) { # name of a hcl palette
      myCols <- hcl.colors(N, pal)
      m <- pal
    } else {   # already a color vector
      myCols <- colorRampPalette(pal)(N)
      m <- sprintf("Custom palette (spread to %d values)", N)
    }
  } else { # palette function
    myCols <- pal(N)
    m <- as.character(substitute(pal))
  }
  barplot(rep(1,N), main = m, axes = FALSE, col =  myCols)
}

plotPal(rainbow)         # rainbow() is the quintessential qualitative palette
plotPal(cm.colors)       # cm.colors() is sequential: cyan-white-magenta
plotPal(terrain.colors)  # a diverging palette: map elevation colors

# A lot of thought has gone into the  construction of palettes in the
# colorspace:: package, now available as hcl.colors in R (hcl: hue, chroma,
# luminance), which is a better way to specify perceptually equidistant colors.
# It's worthwhile to read through the help-file for more information -
?hcl.colors
# ... and to sample the 110 available palettes:
hcl.pals()

# Examples:
plotPal("Viridis")
plotPal("Pastel 1")
plotPal("Spectral")
plotPal("OrRd")
plotPal("Cold")
plotPal("Mint")

# If you are curious, you can execute the code below, then select the line that
# draws the plot and press ctrl+enter repeatedly to step through the entire
# set of palettes.

# i <- 0
# i<-i+1; plotPal(hcl.pals()[i])

# Worthy of special mention are the color-Brewer palettes, in particular for
# their accessibility considerations. Do consider that not all of us view colors
# in the same way!
if (! requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
?RColorBrewer
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)

plotPal(RColorBrewer::brewer.pal(11, "PuOr"), N = 11)


# ===   06.3.2  Custom palettes
#
# Bespoke palettes are easily created with colorRampPalette(). The function
# returns a palette, i.e. a function (not a vector of colors). You assign the
# function to a variable (your palette), and you call it with an argument of how
# many color values you want it to return.


myPal <- colorRampPalette(c("#000000", "#FFFFFF", "#FF0000"))
myPal(10)
plotPal(myPal)

# The parameter "bias" controls the spacing of colors at the end of the palette:
plotPal(colorRampPalette(c("#000000", "#FFFFFF", "#FF0000"), bias = 0.5))
plotPal(colorRampPalette(c("#000000", "#FFFFFF", "#FF0000"), bias = 1.0))
plotPal(colorRampPalette(c("#000000", "#FFFFFF", "#FF0000"), bias = 2.0))

# But a similar effect can be obtained by duplicating one or more colours to
# spread the range ...
plotPal(colorRampPalette(c("#000000", "#FFFFFF", "#FF0000", "#FF0000")))

# ... and this can also be used to spread the middle
plotPal(colorRampPalette(c("#000000", "#FFFFFF", "#FFFFFF", "#FF0000")))

# nicer pastels
myPastels <- colorRampPalette(c("#FFFFFF",
                                "#FF9AA2",
                                "#FFB7B2",
                                "#FFDAC1",
                                "#E2F0CB",
                                "#B5EAD7",
                                "#C7CEEA",
                                "#FFFFFF"))
plotPal(myPastels)

# outliers
myOutliers <- colorRampPalette(c("#66666A",
                                 "#9999A6",
                                 "#CCCCCA",
                                 "#EEEEF9",
                                 "#FFBBA8",
                                 "#F30000"), bias = 0.5)
plotPal(myOutliers)

# Other many sources on the Web that help to generate matching colors upon which
# to build palette. One such site with a very large selection of
# community-contributed palettes is:  https://color.adobe.com/  ... and there is
# also an option to extract palettes from user-supplied images.


# ===   06.3.3  Transparency: The Alpha Channel

# R colours are actually specified as quartets: the fourth value
# the "Alpha channel" defines the transparency. Setting this to
# values other than "FF" (the default) can be useful for very
# crowded plots, or for creating overlays.


x <- SC$xpr[, "t.10"] # towards the beginning of the cycle
y <- SC$xpr[, "t.50"] # towards the end of the cycle

# compare:
plot(x,y, pch = 19, col = "#EE3A8C")
plot(x,y, pch = 19, col = "#EE3A8C12") # Alpha at ~ 10%

# or with multiple overlays of varying size using points() ...
plot(  x,y, pch = 16, cex = 1,   col = "#AA330009")
points(x,y, pch = 19, cex = 2,   col = "#44558803")
points(x,y, pch = 20, cex = 0.5, col = "#EE3A8C08")

# A similar behaviour can be obtained from "density adapted colors with the
# densCols() function
plot (x, y, col = densCols(x, y), pch = 19, cex = 1.5)


# ... or with the smootScatter function
smoothScatter(x, y)

smoothScatter(x, y,
              colramp = colorRampPalette(c("#66666A",
                                           "#9999A6",
                                           "#CCCCCA",
                                           "#EEEEF9",
                                           "#FFBBA8",
                                           "#F30000"), bias = 0.5),
              nrpoints = 0,
              xlim = c(-0.4, 0.3),
              ylim = c(-0.4, 0.3),
              xlab = "Expression (t = 10 min.)",
              ylab = "Expression (t = 50 min.)")
abline(h = 0, col = "#FFFFFF", lwd = 0.5)
abline(v = 0, col = "#FFFFFF", lwd = 0.5)

# ToDo: expand to "real" density plots (cf. )

# ==   06.4  abline(), lines()  and segments()  ================================

# lines() draws an arbitrary line from the supplied points

N <- 300
x <- seq(0.8 * pi, 40 * pi, length.out = N)
y <- sin(x)/x
plot(x, y, type = "n") # only set up the axes
abline(h = 0)
lines(x, y, lwd = 4, col = "#FFFFFF")  # "erase" parts of the axis
lines(x, y, col = "#CC0000")           # draw the line we want


# Varying the color of a plotted line can't be done in base r. Use
# ggplot2:: or plotrix::color.scale.lines() instead

if (! requireNamespace("plotrix", quietly = TRUE)) {
    install.packages("plotrix")
}

# plotrix::color.scale.lines() can change color as well as width of a line
# to give you more ways to visually display features of your data

N <- 300
x <- seq(0.8 * pi, 40 * pi, length.out = N)
y <- sin(x)/x
myCol <- myPastels(N)                   # colors
myLwd <- seq(0.5, 10, length.out = N)    # line width

plot(x, y, type = "n") # only set up the axes
abline(h = 0, lwd = 2, col = "#CCCCCC")

plotrix::color.scale.lines(x, y, col = myCol, lwd = myLwd)

# =    07  AXES  ===============================================================

# For Details, see:
?plot.default

n <- 1000
x <- rnorm(n)
y <- x^3 * 0.25 + rnorm(n, sd=0.75)

plot(x,y)  # Default

# Axes
plot(x,y, xlim=c(-4, 4)) # fixed limits
plot(x,y, xlim=c(-4, 4), ylim=c(10, -10)) # reverse is possible
plot(x,y, log="xy") # log axes

# The axis parameters in the default plot are limited.
# If you want more control, suppress the printing of an axis
# in the plot and use the axis() function instead.
?axis


# Axis-labels and title are straightforward parameters of plot
plot(x,y, xlab="rnorm(n)",
          ylab="x^3 * 0.25 + rnorm(n, sd=0.75)",
          cex.main=1.3,
          main="Sample\nPlot",
          cex.sub=0.75,
          col.sub="grey",
          sub="Scatterplot of noisy 3d-degree polynomial"
          )

# Add gridlines
?grid
grid()


# =    08  LEGENDS  ============================================================


# ==   08.1  basic legends  ====================================================
#
# ToDo: basic legend syntax

# ==   08.2  Color bars  =======================================================

# Legends for a plot with a sequential or divergent spectrum are a special case, since we need to map a continuum of colors to a scale

# Example: plotting expression levels in two cell-cycle phases against each other, and applying a divergent color spectrum that displays the model correlations: genes with low correlations have expression profiles that can not be well approximated with a periodic function.


tRep <- c("t.15", "t.20", "t.25")   # Replication-phase time points
tDup <- c("t.40", "t.45", "t.50")   # Duplication-phase time points
x <- rowMeans(SC$xpr[ , tRep])
y <- rowMeans(SC$xpr[ , tDup])
plot(x, y)

myCors <- SC$mdl$cor
hist(myCors, breaks = 40)

# The values are of course continuous, but we need to map them to a discrete set of colors. The strategy is as follows:
#  - scale the values we wish to map into the interval [0, 1]
#  - multiply with the desired number of intervals minus one, round, and add one
#  - these integers can  be used as indices into a vector of colors, ...
#  - ... which gives the colors for the plot

# A divergent palette
corPal <- colorRampPalette(c("#90005DAA", "#D44F99AA", "#F4C2B8AA",
                             "#F8DEB9AA", "#EEFAE9AA", "#D7F1C8AA",
                             "#BBE9C4AA", "#7CB8A6AA"), alpha = TRUE)
hist(myCors,
     breaks = seq(min(myCors), max(myCors), length.out = 21),
     col = corPal(20))

N <- 10  # We'll map into N  intervals


#$$$
# The usual way to scale a vector x to the [0,1] interval
# is:  x-min / max-min
idx <- (myCors - min(myCors)) / (max(myCors) - min(myCors))
# Now: expand this to N-1, round, and add 1
idx <- round(idx * (N-1)) + 1
# confirm that these numbers are now suitable as indices into an 1:N vector
table(idx)
# create the color vector
( corCols <- corPal(N) )
# and make a vector of colors, one for each gene
myCols <- corCols[idx]

# then we plot ...
plot(x, y, pch = 16, col = myCols, cex = 1.2)

# there is a lot of overlap - let's re-plot in order of decreasing correlations
# so the low correlations (reds) get plotted last. Ordering this way will give a
# "cleaner" appearance for divergent plots.
ord <- order(myCors, decreasing = TRUE)
plot(x[ord], y[ord], pch = 16, col = myCols[ord], cex = 1.2)

#: Titles etc.
basePlot <- function() {
  m <- sprintf("Expression changes between phases (%d genes)",
               length(x))
  plot(x[ord], y[ord],
       main = m,
       xlab = "mean Expression in Replication phase",
       ylab = "mean Expression in Duplication phase",
       sub = paste("Coloured by Correlation of",
                   "Expression Profile (xpr) to",
                   "Fitted Periodic Function (mdl)."),
       cex.main = 0.9,
       cex.lab = 0.7,
       cex.sub = 0.7,
       pch = 16,
       col = myCols[ord],
       cex = 1.2)
  abline(h = 0, col = "#0088FF44")       # horizontal zero
  abline(v = 0, col = "#0088FF44")       # vertical zero
}

basePlot()

# to print a color-bar legend, we use the plotrix package

if (! requireNamespace("plotrix", quietly = TRUE)) {
  install.packages("plotrix")
}

corLabels <- sprintf("%4.2f",           # craft the labels to plot
                     seq(min(myCors),
                         max(myCors),
                         length.out = length(corCols)))

oPar <- par(mar=c(5,4,4,7))  # more space on the margin
par("xpd" = FALSE)
basePlot()

# we need to figure out where to put the legend. par("usr") gives us the coordinates of the last plot window. We'll space the legend between 1.05 and 1.1 times the plot width, and centred along the y-axis with a height of 0.8 times the plot height.
dx <- par("usr")[2] - par("usr")[1]
dy <- par("usr")[4] - par("usr")[3]
xl <- par("usr")[2] + (dx * 0.05)    # left
xr <- xl + (0.05 * dx)               # right
yb <- par("usr")[4] - (0.9 * dy)     # bottom
yt <- yb + (0.8 * dy)                # top


plotrix::color.legend(xl, yb, xr, yt,
                      corLabels,
                      corCols,
                      align = "rb",       # labels to the right of the bar
                      cex = 0.8,
                      gradient = "y")     # vertical gradient

# Almost perfect. Except we need a title for the color bar.
# in this case the title is "r:" ... but the "r" should be in italics.
par("xpd" = TRUE)  # enable plotting into the margin
text(xr + (0.01 * dx), yt + 0.015 * dy,
     expression(italic(r)["xpr,mdl"]),
     adj = c(0, 0),
     cex = 0.9)

par(oPar)

# ToDo - this is a common plot prototype - cast it into a function

# =    09  LAYOUT  =============================================================

# par, lattice, constant aspect ratio


# Most parameters of the plot window can be set via
# the functions plot(), hist() etc., but some need to
# be set via the par() function. Calling par() without
# arguments lists the current state of the plotting
# parameters. Calling it with arguments, returns the
# old parameters and sets new parameters. Thus setting
# new parameters and saving the old ones can be done
# in one step. The parameters that have to be set via
# par include:

# -  multiple plots in one window (mfrow, mfcol, mfg)
# - margin layout (mai, mar mex, oma, omd, omi)
# - controlling position and size of a plot in the figure (fig, plt, ps, pty)
# - see ?par for details.

n <- 1000
x <- rnorm(n)
y <- x^3 * 0.25 + rnorm(n, sd=0.75)

# set window background and plotting axes via par
opar <- par(bg="steelblue", fg="lightyellow")
# set axis lables and titles via plot parameters
plot(x,y, col.axis="lightyellow", col.lab="lightyellow")
par(opar)  # reset to old values

plot(x,y) # confirm reset

# aspect! ...
#

# =    10  TEXT  ===============================================================


# ToDo: sprintf()

# Plotting arbitrary text
# use the text() function to plot characters and strings to coordinates
?text

# Example: add labels to the symbols
# first set: plain symbols (1 to 20)
text(x1-0.4, y1, paste(1:20), cex=0.75)
# symbols with separate background (21 to 25)
text(x2-0.4, y2, paste(21:25), cex=0.75)
# third set: special characters, change font for clarity
text(x3-0.4, y3, extra, col="slateblue", cex=0.75, vfont=c("serif", "plain"))


# ToDo: writing into the margin ...
#   par("mar")
#   par("xpd")
#   par("srt")
#
# ToDo: plotmath()   # demo(plotmath)
# ToDo: expression()



# =    11  DRAWING ON PLOTS  ===================================================


# To plot lines on top of plots there are several options:

# segments()  # for lines with known endpoints
# lines()     # joining x,y coordinate lists with lines
# arrows()    # ... but to get a filled arrow use polygon()
# curve()     # drawing a curve based on a given expression
# abline()    # drawing horizontal or vertical lines or lines with
#             #   given slope and intercept ()


# Example: dividing a plot into 60° regions, centred on a point.
# A general approach to "lines" on a plot is provided by segments().
# However in this special case one can use abline().
# We have to take care though that the aspect ratio for the
# plot is exactly 1 - otherwise our angles are not right.
# Therefore we need to set the asp parameter for plots.

# For a general sketch
#  - we plot the frame a bit larger, don't draw axes
#  - draw the ablines
#  - draw two arrows to symbolize the coordinate axes

p <- c(4, 2)
plot(p[1], p[2],
     xlim=c(-0.5,10.5),
     ylim=c(-0.5,10.5),
     xlab="", ylab="",
     axes=FALSE,
     asp=1.0)
abline(h=p[2], lty=2)  # horizontal
abline(p[2] - (p[1]*tan(pi/3)),  tan(pi/3), lty=2)  # intercept, slope
abline(p[2] + (p[1]*tan(pi/3)), -tan(pi/3), lty=2)  # intercept, slope
arrows(0, 0, 10, 0, length=0.1)   # length of arrow
arrows(0, 0, 0, 10, length=0.1)


# curves()
# rect()
# polygon()

# example: plot the area under the normal distribution as a polygon, and overlay
# a histogram
set.seed(12357)
x <- rnorm(50)
sig <- 1.0
myBreaks <- seq(-3 * sig, 3 * sig, by = 0.5 * sig)

# first, draw a histogram but don't show the bars
hist(x, breaks = myBreaks, freq = FALSE, col = "#FFFFFF", border = "#FFFFFF")

# next, plot the curve as a polygon() into the existing plot area
dn <- seq(-3, 3, length.out = 100)
polygon(dn, dnorm(dn), col = "#00FF8822", border = "#FFFFFF")

# finally, replot the histogram. Setting "add = TRUE" plots it over the
# previous histogram
hist(x, breaks = myBreaks, freq = FALSE,
     col = "#0000FF22", border = "#00000033",   # note: transparent colours
     add = TRUE)

#
# More: see the Index of functions for the graphics package


# =    12  IMAGES  =============================================================

# Todo


# =    13  CONTOUR LINES  ======================================================

# ToDo


# =    14  3D PLOTS  ===========================================================

# ToDo


# =    15  GRAPHS AND NETWORKS  ================================================

# ToDo


# =    16  OTHER GRPAHICS PACKAGES  ============================================

# Packages in the standard distribution ...
#
#   graphics::
#   grid::
#   lattice::

# Packages that can be downloaded from CRAN
# ... use with install.packages("package"), then
#              library("package")

#   hexbin::
#   ggplot2::

# Packages that can be downloaded  from BioConductor
#   prada:
# if (! requireNamespace("prada", quietly=TRUE)) {
#     if (! requireNamespace("BiocManager", quietly=TRUE)) {
#       install.packages("BiocManager")
#     }
#     BiocManager::install("prada")
# }

# =    17  INTERACTIVE PLOTS  ==================================================

# ToDo

# ==   17.1  locator()  ========================================================
# ToDo

# ==   17.2  plotly::  =========================================================
# ToDo




# [End]

