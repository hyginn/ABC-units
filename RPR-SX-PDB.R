# tocID <- "RPR-SX-PDB.R"
#
# ---------------------------------------------------------------------------- #
#  PATIENCE  ...                                                               #
#    Do not yet work wih this code. Updates in progress. Thank you.            #
#    boris.steipe@utoronto.ca                                                  #
# ---------------------------------------------------------------------------- #
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-SX-PDB unit.
#
# Version:  1.2
#
# Date:     2017  10  -  2019  11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    Maintenance
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
#           1.0    First live version, completely refactores 2016 code
#                     with remarkable speed gains. Added section on x, y, z
#                     (density) plots.
#           0.1    First code copied from 2016 material.
#
# TODO:
#          Confirm that SS residue numbers are indices
#          Set task seed from student number
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
#TOC>   Section  Title                                      Line
#TOC> ----------------------------------------------------------
#TOC>   1        Introduction to the bio3D package            62
#TOC>   2        A Ramachandran plot                         153
#TOC>   3        Density plots                               229
#TOC>   3.1        Density-based colours                     243
#TOC>   3.2        Plotting with smoothScatter()             262
#TOC>   3.3        Plotting hexbins                          277
#TOC>   3.4        Plotting density contours                 305
#TOC>   3.4.1          ... as overlay on a coloured grid     338
#TOC>   3.4.2          ... as filled countour                355
#TOC>   3.4.3          ... as a perspective plot             386
#TOC>   4        cis-peptide bonds                           404
#TOC>   5        H-bond lengths                              419
#TOC>
#TOC> ==========================================================================


# In this example of protein structure interpretation, we ...
#   - load the library "bio3D" which supports work with
#     protein structure files,
#   - explore some elementary functions of the library
#   - explore plotting of density values with scatterplots
#   - assemble a script to see whether H-bond lengths systematically differ
#     between alpha-helical and beta-sheet structures.


# =    1  Introduction to the bio3D package  ===================================


if (! requireNamespace("bio3d", quietly = TRUE)) {
  install.packages("bio3d")
}
# Package information:
#  library(help = bio3d)       # basic information
#  browseVignettes("bio3d")    # available vignettes
#  data(package = "bio3d")     # available datasets


# bio3d can load molecules directly from the PDB servers, you don't _have_ to
# store them locally, but you could.
#
# But before you _load_ a file, have a look what such a file contains. I have
# packaged the 1BM8 file into the project:
file.show("./data/1BM8.pdb")

# Have a look at the header section, the atom records, the coordinate records
# etc.
#
# Task: Answer the following questions:
#
#          What is the resolution of the structure?
#          Is the first residue in the SEQRES section also the first residue
#              with an ATOM record?
#          How many water molecules does the structure contain?
#          Can you explain REMARK 525 regarding HOH 459 and HOH 473?
#          Are all atoms of the N-terminal residue present?
#          Are all atoms of the C-terminal residue present?

apses <- bio3d::read.pdb("1bm8")  # load a molecule directly from the PDB via
                                  # the Internet.

# check what we have:
apses

# what is this actually?
str(apses)

# bio3d's pdb objects are simple lists. Great! You know lists!

# You see that there is a list element called $atom which is a data frame in
# which the columns arevectors of the same length - namely the number of atoms
# in the structure file. And there is a matrix of (x, y, z) triplets called xyz.
# And there is a vector that holds sequence, and two tables called helix and
# sheet. Let's pull out a few values to confirm how selection and subsetting
# works here:

# selection by atom ...
i <- 5
apses$atom[i,]
apses$atom[i, c("x", "y", "z")]   # here we are selecting with column names!
apses$xyz[c(i*3-2, i*3-1, i*3)]   # here we are selcting with row numbers

# all atoms of a residue ...
i <- 48
apses$atom[apses$atom[,"resno"] == i, ]

# sequence of the first ten residues
apses$seqres[1:10]  # the "A"s here identify chain "A"

# Can we convert this to one letter code?
bio3d::aa321(apses$seqres[1:10])

# Lets get the implicit sequence:
bio3d::aa321((apses$atom$resid[apses$calpha])[1:10])
# Do you understand this code?

# Do explicit and implicit sequence have the same length?
length(apses$seqres)
length(apses$atom$resid[apses$calpha])

# Are the sequences the same?
sum(apses$seqres == apses$atom$resid[apses$calpha])

# get a list of all CA atoms of arginine residues
sel <- apses$atom$resid == "ARG" & apses$atom$elety == "CA"
apses$atom[sel, c("eleno", "elety", "resid", "chain", "resno", "insert")]

# The introduction to bio3d tutorial at
#   http://thegrantlab.org/bio3d/tutorials/structure-analysis
# has the following example:
bio3d::plot.bio3d(apses$atom$b[apses$calpha],
                  sse=apses,
                  typ="l",
                  ylab="B-factor")

# Good for now. Let's do some real work.

# =    2  A Ramachandran plot  =================================================

# Calculate a Ramachandran plot for the structure. The torsion.pdb() function
# calculates all dihedral angles for backbone and sidechain bonds, NA where
# the bond does not exist in an amino acid.
tor <- bio3d::torsion.pdb(apses)
plot(tor$phi, tor$psi,
     xlim = c(-180, 180), ylim = c(-180, 180),
     main = "Ramachandran plot for 1BM8",
     xlab = expression(phi),
     ylab = expression(psi))
abline(h = 0, lwd = 0.5, col = "#00000044")
abline(v = 0, lwd = 0.5, col = "#00000044")
# As you can see, there are a number of points in the upper-right
# quadrant of the plot. This combination of phi-psi angles defines
# the conformation of a left-handed alpha helix and is generally
# only observed for glycine residues. Let's replot the data, but
# colour the points for glycine residues differently. First, we
# get a vector of glycine residue indices in the structure:

mySeq <- bio3d::pdbseq(apses)

# Explore the result object and extract the indices of GLY resiues.
              mySeq == "G"
        which(mySeq == "G")
iGly <- which(mySeq == "G")

# Now plot all non-gly residues.
# Remember: negative indices exclude items from a vector
plot(tor$phi[-iGly], tor$psi[-iGly],
     xlim=c(-180,180), ylim=c(-180,180),
     main = "Ramachandran plot for 1BM8",
     xlab = expression(phi),
     ylab = expression(psi))
abline(h = 0, lwd = 0.5, col = "#00000044")
abline(v = 0, lwd = 0.5, col = "#00000044")

# Now plot GLY only, but with green dots:
points(tor$phi[iGly], tor$psi[iGly], pch=21, cex=0.9, bg="#00CC00")

# As you see, four residues in the upper-right quadrant are
# not glycine. But what residues are these? Is there an
# error in our script? Let's get their coordinate records:

# subset CA records
CA <- apses$atom[apses$calpha, c("eleno", "elety", "resid", "chain", "resno")]

# get index of outliers
iOutliers <- which(tor$phi > 30 & tor$phi < 90 &
                     tor$psi > 0 & tor$psi < 90)

# cbind records together
(dat <- cbind(CA[iOutliers, ],
              phi = tor$phi[iOutliers],
              psi = tor$psi[iOutliers]))


# remove the glycine ...
(dat <- dat[dat$resid != "GLY", ])


# let's add the residue numbers to the plot with the text function:
for (i in 1:nrow(dat)) {
  points(dat$phi[i], dat$psi[i], pch=21, cex=0.9, bg="#CC0000")
  text(dat$phi[i],
       dat$psi[i],
       labels = sprintf("%s%d", bio3d::aa321(dat$resid[i]), dat$resno[i]),
       pos = 4,
       offset = 0.4,
       cex = 0.7)
}

# You can check the residues in Chimera. Is there anything unusual about these
# residues?


# =    3  Density plots  =======================================================

# Such x, y scatter-plots of data that is sampled from a distribution can tell
# us a lot about what shapes the distribution. The distribution is governed by
# the free energy of the phi-psi landscape in folded proteins, since folded
# proteins generally minimize the free energy of their conformations. We observe
# empirically, from comparing frequency statistics and mutation experiments,
# that this generall follows a Boltzmann distribution, where the free energy
# changes we observe in experments that change one conformation into another are
# proportional to the log-ratio of the number of times we observe each
# observation in the protein structure database (after correcting for
# observation bias). The proper way to visualize such 2D landscapes is with
# contour plots.

# ==   3.1  Density-based colours  =============================================

# A first approximation to  scatterplots that visualize the density of the
# underlying distribution is colouring via the densCols() function.
?densCols
iNA <- c(which(is.na(tor$phi)), which(is.na(tor$psi)))
phi <- tor$phi[-iNA]
psi <- tor$psi[-iNA]
plot (phi, psi,
      xlim = c(-180, 180), ylim = c(-180, 180),
      col=densCols(phi,psi),
      pch=20, cex=2,
      main = "Ramachandran plot for 1BM8",
      xlab = expression(phi),
      ylab = expression(psi))
abline(h = 0, lwd = 0.5, col = "#00000044")
abline(v = 0, lwd = 0.5, col = "#00000044")


# ==   3.2  Plotting with smoothScatter()  =====================================

# A better way, that spreads out the underlying density is smoothScatter()
smoothScatter(phi,psi)
smoothScatter(phi, psi,
              xlim = c(-180, 180), ylim = c(-180, 180),
              col = "#0033BB33",
              pch = 3, cex = 0.6,
              main = "Ramachandran plot for 1BM8",
              xlab = expression(phi),
              ylab = expression(psi))
abline(h = 0, lwd = 0.5, col = "#00000044")
abline(v = 0, lwd = 0.5, col = "#00000044")


# ==   3.3  Plotting hexbins  ==================================================

# If we wish to approximate values in a histogram-like fashion, we can use
# hexbin()
if (! requireNamespace("hexbin", quietly = TRUE)) {
  install.packages("hexbin")
}
# Package information:
#  library(help = hexbin)       # basic information
#  browseVignettes("hexbin")    # available vignettes
#  data(package = "hexbin")     # available datasets


myColorRamp <- colorRampPalette(c("#EEEEEE",
                                  "#3399CC",
                                  "#2266DD"))
hexbin::gplot.hexbin(hexbin::hexbin(phi, psi, xbins = 10),
                     colramp = myColorRamp,
                     main = "phi-psi Density Bins for 1BM8",
                     xlab = expression(phi),
                     ylab = expression(psi))

# Note:
# Had we loaded hexbin with library(hexbin), the plot function would have
# been dispatched by the plot() generic, and we could simply have written:
#   plot(hexbin(phi, psi, xbins = 10), ...


# ==   3.4  Plotting density contours  =========================================


# The best way to handle such data is provided by the function contour():
?contour

# Contour plots are not produced along the haphazardly sampled values of a data
# set, but on a regular grid. This means, we need to convert observed values
# into estimated densities. Density estimation is an important topic for
# exploratory data analysis, base R has the density() function for 1D
# distributions. But for 2D data like or phi-psi plots, we need a function from
# the MASS package: kde2d()

if (! requireNamespace("MASS", quietly = TRUE)) {
  install.packages("MASS")
}
# Package information:
#  library(help = MASS)       # basic information
#  browseVignettes("MASS")    # available vignettes
#  data(package = "MASS")     # available datasets

?MASS::kde2d
dPhiPsi <-MASS::kde2d(phi, psi,
                      n = 60,
                      lims = c(-180, 180, -180, 180))

str(dPhiPsi)
# This is a list, with gridpoints in x and y, and the estimated densities in z.

# Generic plot with default parameters
contour(dPhiPsi)


# ===   3.4.1  ... as overlay on a coloured grid

image(dPhiPsi,
      col = myColorRamp(100),
      main = "Ramachandran plot for 1BM8",
      xlab = expression(phi),
      ylab = expression(psi))
contour(dPhiPsi, col = "royalblue",
        add = TRUE,
        method = "edge",
        nlevels = 10,
        lty = 2)
points(phi, psi, col = "#00338866", pch = 3, cex = 0.7)
abline(h = 0, lwd = 0.5, col = "#00000044")
abline(v = 0, lwd = 0.5, col = "#00000044")


# ===   3.4.2  ... as filled countour

filled.contour(dPhiPsi,
               xlim = c(-180, 180), ylim = c(-180, 180),
               nlevels = 10,
               color.palette = myColorRamp,
               main = "Ramachandran plot for 1BM8",
               xlab = expression(phi),
               ylab = expression(psi))

# Note: we can pass additional plotting and overlay commands to the counter plot
# in a block of expressions passed via the plot.axes parameter:

filled.contour(dPhiPsi,
               xlim = c(-180, 180), ylim = c(-180, 180),
               nlevels = 10,
               color.palette = myColorRamp,
               main = "Ramachandran plot for 1BM8",
               xlab = expression(phi),
               ylab = expression(psi),
               plot.axes = {
                  contour(dPhiPsi, col = "#00000044",
                          add = TRUE,
                          method = "edge",
                          nlevels = 10,
                          lty = 2)
                    points(phi, psi, col = "#00338866", pch = 3, cex = 0.7)
                    abline(h = 0, lwd = 0.5, col = "#00000044")
                    abline(v = 0, lwd = 0.5, col = "#00000044")
                  })

# ===   3.4.3  ... as a perspective plot

persp(dPhiPsi,
      xlab = "phi",
      ylab = "psi",
      zlab = "Density")


persp(dPhiPsi,
      theta = 40,
      phi = 10,
      col = "#99AACC",
      xlab = "phi",
      ylab = "psi",
      zlab = "Density")



# =    4  cis-peptide bonds  ===================================================

# Are there any cis-peptide bonds in the structure?
tor$omega
#... gives us a quick answer. But wait - what values do we
# expect? And why are the values so different, ranging from -180° to 180°?
# Consider this plot: what am I doing here and why?
om <- c(360 + tor$omega[tor$omega < 0],
        tor$omega[tor$omega >= 0])
hist(om, xlim=c(0,360))
abline(v=180, col="red")

# Note: a cis-peptide bond will have an omega torsion angle around 0°


# =    5  H-bond lengths  ======================================================

# Let's do something a little less trivial and compare
# backbone H-bond lengths between helices and strands.

# Secondary structure is defined in the bio3d object's ...$helix and ...$strand
# list elements.
str(apses)

# We need to
# - collect all N atoms in alpha helices resp.
#      beta strands;
# - collect all O atoms in alpha helices resp.
#      beta strands;
# - fetch the atom coordinates;
# - calculate all N, O distances;
# - filter them for distances that indicate H-bonds; and,
# - plot the results.

# Secondary structure can either be obtained from definitions contained in most
# PDB files, or by running the DSSP algorithm IF(!) you have it installed on
# your machine. See the dssp() function of bio3d for instructions how to obtain
# and install DSSP and STRIDE. This is highly recommended for "real" work with
# structure coordinate files. The 1BM8 file contains secondary structure
# definitions:

apses$helix
apses$sheet


# A function to collect atom indices for particular type of secondary structure

ssSelect <- function(PDB, myChain = "A", ssType, myElety) {
  # Function to retrieve specified atom types from specified secondary
  # structure types in a PDB object.
  # Parameters:
  #    PDB              A bio3D PDB object
  #    myChain     chr  The chain to use. Default: chain "A".
  #    ssType      chr  A vector of keywords "helix" and/or "sheet"
  #    myElety     chr  A vector of $eletype atom types to return
  # Value:         num  Indices of the matching atom rows in PDB$atom

  # Build a vector of $resno numbers
  starts <- numeric()
  ends   <- numeric()
  if ("helix" %in% ssType) {
    sel <- PDB$helix$chain %in% myChain
    starts <- c(starts, PDB$helix$start[sel])
    ends   <- c(ends,   PDB$helix$end[sel])
  }
  if ("sheet" %in% ssType) {
    sel <- PDB$sheet$chain %in% myChain
    starts <- c(starts, PDB$sheet$start[sel])
    ends   <- c(ends,   PDB$sheet$end[sel])
  }
  myResno <- numeric()
  for (i in seq_along(starts)) {
    myResno <- c(myResno, starts[i]:ends[i])
  }

  # get id's from PDB

  x <- bio3d::atom.select(PDB,
                          string = "protein",
                          type = "ATOM",
                          chain = myChain,
                          resno = myResno,
                          elety = myElety)

  return(x$atom)
}

# Example:
ssSelect(apses, ssType = "sheet", myElety = "N")
ssSelect(apses, ssType = "sheet", myElety = "O")

# That looks correct: O atoms should be stored three index position after N: the
# sequence of atoms in a PDB file is usually N, CA, C, O ... followed by the
# side chain coordinates.

# Now to extract the coordinates and calculate distances. Our function needs to
# take the PDB object and two vectors of atom indices as argument, and return a
# vector of pair-distances (actually  dist.xyz() returns a matrix).

pairDist <- function(PDB, a, b) {
  # Calculate pairwise distances between atoms indicated by a and b
  # Parameters:
  #    PDB       A bio3D PDB object
  #    a    int  A vector of atom indexes
  #    b    int  A vector of atom indexes
  # Value:  num  matrix of Euclidian distances between the atoms given in a, b.
  #                There are as many rows as atoms in a, as many columns as
  #                atoms in b.

  dMat <- numeric()
  if (length(a) > 0 && length(b) > 0) {

  A <- PDB$atom[a, c("x", "y", "z")]
  B <- PDB$atom[b, c("x", "y", "z")]
  dMat <- bio3d::dist.xyz(A, B)

  }
  return(dMat)
}

# Let's see if this looks correct. Let's look at all the pairwise distances
# between N and O atoms in both types of secondary structure:

iN <- ssSelect(apses, ssType = c("helix", "sheet"), myElety = "N")
iO <- ssSelect(apses, ssType = c("helix", "sheet"), myElety = "O")
x <- pairDist(apses, iN, iO)
hist(x,
     breaks = 80,
     col = "lavenderblush",
     main = "",
     xlab = "N-O distances (Å)")

# Since we are collecting distance from all secondary structure elements, we
# are just seing a big peak of (meaningless) long-distance separations. We
# need to zoom in on the shorter distances, in which we expect
# hydrogen bonds:
hist(x[x < 4.2],                # restrict to N-O distance less than 4.2 Å long
     breaks=seq(2.0, 4.2, 0.1),
     xlim=c(2.0,4.2),
     col = "lavenderblush",
     main = "",
     xlab = "N-O distances (Å)")

# There is a large peak at about 2.2Å, and another
# large peak above 3.5Å. But these are not typical hydrogen
# bond distances! Rather these are (N,O) pairs in peptide
# bonds, and within residues. That's not good, because these will contaminate
# our statistics.
# We need to exclude all distances between N of a residue
# and O of a preceeding residue, and all (N,O) pairs in the
# same residue. We need a function to filter distances by residue numbers. And while we are filtering, we might as well throw away the non-H bond distances too.

filterHB <- function(PDB, iN, iO, dMat, cutoff = 3.9) {
  # Filters distances between O(i-1) and N(i), and between N(i) and O(i)
  # in a distance matrix where there is one row per N-atom and one
  # column per O atom.
  # Parameters:
  #    PDB            a bio3D PDB object
  #    iN       int   a vector of N atom indexes
  #    iO       int   a vector of O atom indexes
  #    dMat     num   a distance matrix created by pairDist()
  #    cutoff   num   only return distances that are shorter than "cutoff
  # Value:  a distance matrix in which values that do not match the
  #            filter criteria have bee set to NA.

  if (length(iN) > 0 && length(iO) > 0) {

    resN <- PDB$atom$resno[iN]
    resO <- PDB$atom$resno[iO]

    for (i in seq_along(resN)) {        # for all N atoms
      for (j in seq_along(resO)) {      # for all O atoms
        if (dMat[i, j] > cutoff ||      # if: distance > cutoff, or ...
            (resN[i] - 1) == resO[j] || #     distance is N(i)---O(i-1), or ...
            resN[i] == resO[j]) {       #     distance is N(i)---O(i), then:
          dMat[i, j] <- NA              # set this distance to NA.
        }
      }
    }
  }
  return(dMat)
}

# Inspect the result:
hist(filterHB(apses, iN, iO, x),
     breaks=seq(2.0, 4.2, 0.1),
     xlim=c(2.0,4.2),
     col = "paleturquoise",
     main = "",
     xlab = "N-O distances (Å)")


# Finally we can calculate alpha- and beta- structure
# bonds and compare them. In this section we'll explore
# different options for histogram plots.

# H-bonds in helices ...
iN <- ssSelect(apses, ssType = c("helix"), myElety = "N")
iO <- ssSelect(apses, ssType = c("helix"), myElety = "O")
dH <- filterHB(apses, iN, iO, pairDist(apses, iN, iO))
dH <- dH[!is.na(dH)]

# H-bonds in sheets. (We commonly use the letter "E" to symbolize a beta
# strand or sheet, because "E" visually evokes an extended strand with
# protruding sidechains.)
iN <- ssSelect(apses, ssType = c("sheet"), myElety = "N")
iO <- ssSelect(apses, ssType = c("sheet"), myElety = "O")
dE <- filterHB(apses, iN, iO, pairDist(apses, iN, iO))
dE <- dE[!is.na(dE)]

# The plain histogram functions without parameters
# give us white stacks.

hist(dH)

# and ...
hist(dE)

# We can see that the histrograms look different
# but that is better visualized by showing two plots
# in the same window. We use the par() function, for
# more flexible layout, look up the layout() function.
?par
?layout

opar <- par(no.readonly=TRUE)  # store current state
par(mfrow=c(2,1))  # set graphics parameters: 2 rows, one column

# plot two histograms
hist(dH)
hist(dE)


# add colour:
hist(dH, col="#DD0055")
hist(dE, col="#00AA70")



# For better comparison, plot both in the
# same window:

hist(dH, col="#DD0055")
hist(dE, col="#00AA70", add=TRUE)

# ... oops, we dind't reset the graphics parameters. You can either close the
# window, a new window will open with default parameters, or ...
par(opar)      # ... reset the graphics parameters

hist(dH, col="#DD0055")
hist(dE, col="#00AA70", add=TRUE)

# We see that the leftmost column of the sheet bonds hides the helix bonds in
# that column. Not good. But we can make the colours transparent! We just need to
# add a fourth set of two hexadecimal-numbers to the #RRGGBB triplet. Lets use
# 2/3 transparent, in hexadecimal, 1/3 of 256 is x55 - i.e. an RGB triplet
# specied as #RRGGBB55 is only 33% opaque:

hist(dH, col="#DD005555")
hist(dE, col="#00AA7055", add=TRUE)

# To finalize the plots, let's do two more things: Explicitly define the breaks,
# to make sure they match up - otherwise they would not need to... like in this
# example:

hist(dH, col="#DD005555")
hist(dE[dE < 3], col="#00AA7055", add=TRUE)

# Breaks are a parameter in hist() that can either be a scalar, to define how
# many columns you want, or a vector, that defines the actual breakpoints.
brk=seq(2.4, 4.0, 0.1)

hist(dH, col="#DD005555", breaks=brk)
hist(dE, col="#00AA7055", breaks=brk, add=TRUE)

# The last thing to do is to think about rescaling the plot. You notice that the
# y-axis is scaled in absolute frequency (i.e. counts). That gives us some
# impression of the relative frequency, but it is of course skewed by observing
# relatively more or less of one type of secondary structure in a protein. As
# part of the hist() function we can rescale the values so that the sum over all
# is one: set the prameter freq=FALSE.

hist(dH, col="#DD005555", breaks=brk, freq=FALSE)
hist(dE, col="#00AA7055", breaks=brk, freq=FALSE, add=TRUE)

# Adding labels and legend ...

hH <- hist(dH,
           freq=FALSE,
           breaks=brk,
           col="#DD005550",
           xlab="(N,O) distance (Å)",
           ylab="Density",
           ylim=c(0,4),
           main="Helix and Sheet H-bond lengths")
hE <- hist(dE,
           freq=FALSE,
           breaks=brk,
           col="#00AA7060",
           add=TRUE)

legend("topright",
       c(sprintf("alpha (N = %3d)", sum(hH$counts)),
         sprintf("beta  (N = %3d)", sum(hE$counts))),
       fill = c("#DD005550", "#00AA7060"), bty = 'n',
       border = NA)


# With all the functions we have defined,
# it is easy to try this with a larger protein.
# 3ugj for example is VERY large.

pdb <- bio3d::read.pdb("3ugj")

# helices...
iN <- ssSelect(pdb, ssType = c("helix"), myElety = "N")
iO <- ssSelect(pdb, ssType = c("helix"), myElety = "O")
dH <- filterHB(pdb, iN, iO, pairDist(pdb, iN, iO))
dH <- dH[!is.na(dH)]

# sheets
iN <- ssSelect(pdb, ssType = c("sheet"), myElety = "N")
iO <- ssSelect(pdb, ssType = c("sheet"), myElety = "O")
dE <- filterHB(pdb, iN, iO, pairDist(pdb, iN, iO))
dE <- dE[!is.na(dE)]

# histograms ...
brk=seq(2.4, 4.0, 0.1)

hH <- hist(dH,
           freq=FALSE,
           breaks=brk,
           col="#DD005550",
           xlab="(N,O) distance (Å)",
           ylab="Density",
           ylim=c(0,4),
           main="Helix and Sheet H-bond lengths")
hE <- hist(dE,
           freq=FALSE,
           breaks=brk,
           col="#00AA7060",
           add=TRUE)

legend('topright',
       c(paste("alpha (N = ", sum(hH$counts), ")"),
         paste("beta  (N = ", sum(hE$counts), ")")),
       fill = c("#DD005550", "#00AA7060"), bty = 'n',
       border = NA,
       inset = 0.1)

# It looks more and more that the distribution is indeed different. Our sample
# is large, but derives from a single protein. To do database scale statistics,
# we should look at many more proteins. To give you a sense of how, let's do
# this for just ten proteins, randomly selected from non-homologous,
# high-resolution, single domain structures. I have provided a utility function
# for such a selection (details beyond the scope of this project).

myPDBs <- selectPDBrep(10)
# My selection is "2OVJ", "1HQS", "3BON", "4JZX", "3BQ3", "2IUM", "2C9E",
# "4X1F", "2V3I", "3GE3". Yours will be different.

# We are loading the files online - don't do this if you have bandwidth
# limitations.


dH <- c() # collect all helix H-bonds here
dE <- c() # collect all sheet H-bonds here

for (i in seq_along(myPDBs)) {
  pdb <- bio3d::read.pdb(myPDBs[i])

  # helices...
  iN <- ssSelect(pdb, ssType = c("helix"), myElety = "N")
  iO <- ssSelect(pdb, ssType = c("helix"), myElety = "O")
  x  <- filterHB(pdb, iN, iO, pairDist(pdb, iN, iO))
  dH <- c(dH, x[!is.na(x)])

  # sheets
  iN <- ssSelect(pdb, ssType = c("sheet"), myElety = "N")
  iO <- ssSelect(pdb, ssType = c("sheet"), myElety = "O")
  x  <- filterHB(pdb, iN, iO, pairDist(pdb, iN, iO))
  dE <- c(dE, x[!is.na(x)])
}

# Inspect the results

length(dH)  # 4415 (your numbers are different, but there should be many)
length(dE)  # 262

brk=seq(2.0, 4.0, 0.1)

hH <- hist(dH,
           freq=FALSE,
           breaks=brk,
           col="#DD005550",
           xlab="(N,O) distance (Å)",
           ylab="Density",
           ylim=c(0,4),
           cex.main = 0.8,
           main="Helix and Sheet H-bond lengths (10 representative structures)")
hE <- hist(dE,
           freq=FALSE,
           breaks=brk,
           col="#00AA7060",
           add=TRUE)

legend('topright',
       c(paste("alpha (N = ", sum(hH$counts), ")"),
         paste("beta  (N = ", sum(hE$counts), ")")),
       fill = c("#DD005550", "#00AA7060"), bty = 'n',
       border = NA,
       inset = 0.1)

# Task:
#    Why do you think these distributions are different?
#    At what distance do you think H-bonds have the lowest energy?
#    For alpha-helices? For beta-strands?




# [END]
