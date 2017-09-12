# RPR-SX-PDB.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-SX-PDB unit.
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

# = 1 ___Section___

# In this example of protein structure interpretation, we ...
#   - load the library "bio3D" which supports work with
#     protein structure files,
#   - explore some elementary functions of the library
#   - assemble a script to see whether H-bond lengths systematically differ
#     between alpha-helical and beta-sheet structures.

if(!require(bio3d)) {
  install.packages("bio3d", dependencies=TRUE)
  library(bio3d)
}

lbio3d() # ... lists the newly installed  functions,
# they all have help files associated.
# More information is available in the so-called
# "vignettes" that are distributed with most R packages:
vignette("bio3d_vignettes")

# bio3d can load molecules directly from the PDB servers, you don't _have_ to
# store them locally, but you could.
#
# But before you _load_ a file, have a look what such a file contains. I have packaged the 1BM8 file into the project:
file.show("./assets/1BM8.pdb")

# Have a look at the header section, the atom records, the coordinate records
# etc. Answer the following questions:
#
# What is the resolution of the structure?
# Is the first residue in the SEQRES section also the first residue
#     with an ATOM record?
# How many water molecules does the structure contain?
# Can you explain REMARK 525 regarding HOH 459 and HOH 473?
# Are all atoms of the N-terminal residue present?
# Are all atoms of the C-terminal residue present?

apses <- read.pdb("1bm8")  # load a molecule directly from PDB

# check what we have:
apses

# what is this actually?
str(apses)

# bio3d's pdb objects are simple lists. Great! You know lists!

# You see that there is a table called atom in which the elements are vectors of the same length - namely the number of atoms in the structure file. And there is a matrix of (x, y, z) triplets called xyz. And there is a vector that holds sequence, and two tables called helix and sheet that look a lot like our feature annotation tables - in fact many of the principles to store this strutcure data are similar to how we constructed our protein database. Let's pull out a few values to confirm how selection and subsetting works here:

# selection by atom ...
i <- 5
apses$atom[i,]
apses$atom[i, c("x", "y", "z")]   # here we are selecting with column names!
apses$xyz[c(i*3-2, i*3-1, i*3)]   # here we are selcting with row numbers

# all atoms of a residue ...
i <- "48"	#string!
apses$atom[apses$atom[,"resno"] == i, ]

# sequence of the first ten residues
apses$seqres[1:10]  # the As here identify the chain

# Can we convert this to one letter code?
aa321(apses$seqres[1:10])

# Lets get the implicit sequence:
aa321((apses$atom$resid[apses$calpha])[1:10])  # Do you understand this code?

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
plot.bio3d(apses$atom$b[apses$calpha], sse=apses, typ="l", ylab="B-factor")

# Good for now. Let's do some real work.

# ==============================================================================
#        PART TWO: A Ramachandran plot
# ==============================================================================

# Calculate a Ramachandran plot for the structure
tor <- torsion.pdb(apses)
plot(tor$phi, tor$psi)

# As you can see, there are a number of points in the upper-right
# quadrant of the plot. This combination of phi-psi angles defines
# the conformation of a left-handed alpha helix and is generally
# only observed for glycine residues. Let's replot the data, but
# color the points for glycine residues differently. First, we
# get a vector of glycine residue indices in the structure:

sSeq <- pdbseq(apses)

# Explore the result object and extract the indices of GLY resiues.
sSeq == "G"
which(sSeq == "G")
as.numeric(which(sSeq == "G"))
iGly <- as.numeric(which(sSeq == "G"))

# Now plot all non-gly residues.
# Remember: negative indices exclude items from a vector
plot(tor$phi[-iGly], tor$psi[-iGly], xlim=c(-180,180), ylim=c(-180,180))

# Now plot GLY only, but with green dots:
points(tor$phi[iGly], tor$psi[iGly], pch=21, cex=0.9, bg="#00CC00")

# As you see, four residues in the upper-right quadrant are
# not glycine. But what residues are these? Is there an
# error in our script? Let's get their coordinate records:

iOutliers <- which(tor$phi > 30 & tor$phi < 90 &
                     tor$psi > 0 & tor$psi < 90)
CA <- apses$atom[apses$calpha, c("eleno", "elety", "resid", "chain", "resno")]
dat <- cbind(CA[iOutliers, ], phi = tor$phi[iOutliers], psi = tor$psi[iOutliers])
dat

# remove the glycine from our table ...
dat <- dat[dat$resid != "GLY", ]
dat

# let's add the residue numbers to the plot with the text function:
for (i in 1:nrow(dat)) {
  points(dat$phi[i], dat$psi[i], pch=21, cex=0.9, bg="#CC0000")
  text(dat$phi[i],
       dat$psi[i],
       labels = sprintf("%s%d", aa321(dat$resid[i]), dat$resno[i]),
       pos = 4,
       offset = 0.4,
       cex = 0.7)
}

# You can check the residues in Chimera. Is there anything unusual?

# Are there any cis-peptide bonds in the structure?
tor$omega
#... gives us a quick answer. But wait - what values do we
# expect? And why are the values so different?
# Consider this plot: what am I doing here and why?
x <- tor$omega[tor$omega > 0]
x <- c(x, 360 + tor$omega[tor$omega < 0])
hist(x, xlim=c(90,270))
abline(v=180, col="red")



# ==============================================================================
#        PART THREE: H-bond lengths
# ==============================================================================

# Let's do something a little less trivial and compare
# backbone H-bond lengths between helices and strands.
#
# Secondary structure is defined in the list components ...$helix
# and ...$strand.

# We need to
# - collect all residue indices for alpha helices resp.
#      beta strands;
# - fetch the atom coordinates;
# - calculate all N, O distances using dist.xyz();
# - filter them for distances that indicate H-bonds; and,
# - plot the results.

# Secondary structure can either be obtained from
# definitions contained in many PDB files, or by running
# the DSSP algorithm IF(!) you have it installed on your
# machine. The 1BM8 file contains definitions

apses$helix
apses$sheet


# (1): collect the residue numbers
# between the segment boundaries.

H <- numeric() # This will contain the helix residue numbers
for (i in 1:length(apses$helix)) {
  H <- c(H, apses$helix$start[i]:apses$helix$end[i])
}

# Doing the same for the sheet residue numbers requires
# very similar code. Rather than retype the code, it is
# better to write a function that can do both.

getSecondary <- function(sec) {
  iRes <- c()
  for (i in 1:length(sec$start)) {
    iRes <- c(iRes, sec$start[i]:sec$end[i])
  }
  return(iRes)
}




# Compare ...
H
getSecondary(apses$helix)

# ... and use for strands

E <- getSecondary(apses$sheet)


# Now here's a problem: these numbers refer to the
# residue numbers as defined in the atom records. They
# can't be used directly to address e.g. the first, second
# third residue etc. since the first residue has the
# residue number 4...
apses$atom[1,]

# Therefore we need to
# 1: convert the numbers to strings;
# 2: subset the atom table for rows contain these strings.
#
# Essentially, we don't treat the "residue numbers" as numbers,
# but as labels. That's fine, as long as they are unique.

# (2): fetch coordinates of N and O atoms
# for residues in alpha- and beta- conformation.

# For one residue, the procedure goes as follows:

res <- H[17] # pick an arbitrary alpha-helix residue to illustrate
res
res <- as.character(res)
res

# all atom rows for this residue
apses$atom[apses$atom[,"resno"] == res, ]

# add condition: row with "N" atom only
apses$atom[apses$atom[,"resno"] == res &
             apses$atom[,"elety"] == "N", ]

# add column selection: "x", "y", "z"
apses$atom[apses$atom[,"resno"] == res &
             apses$atom[,"elety"] == "N",
           c("x", "y", "z")]

# convert to numbers
as.numeric (
  apses$atom[apses$atom[,"resno"] == res &
               apses$atom[,"elety"] == "N",
             c("x", "y", "z")]
)

# Now we need to add this into a matrix as we iterate over the desired residues.
# We need to execute the procedure four times for alpha and beta Ns and Os
# respectively. Rather than duplicate the code four times over, we write a
# function. Never duplicate code, because if you need to make changes it is too
# easy to forget making the change consistently in all copies.


getAtom <- function(PDB, r, AT) {
  # Function to iterate over residue number strings and
  # return a matrix of x, y, z triplets for each atom
  # of a requested type.
  mat <- c() 	#initialize as empty matrix
  for (i in 1:length(r)) {
    res <- as.character(r[i])
    v <- as.numeric (
      PDB$atom[PDB$atom[,"resno"] == res &
                 PDB$atom[,"elety"] == AT,
               c("x", "y", "z")]
    )
    mat <- rbind(mat, v)
  }
  return(mat)
}

# Now run the functions with the parameters we need...
H.xyz.N <- getAtom(apses, H, "N")  # backbone N atoms in helix
H.xyz.O <- getAtom(apses, H, "O")  # backbone O atoms in helix
E.xyz.N <- getAtom(apses, E, "N")  # backbone N atoms in strand
E.xyz.O <- getAtom(apses, E, "O")  # backbone O atoms in strand


# (3): calculate distances between N and O atoms to find H-bonds.

# We spent most of our effort so far just preparing our raw data for analysis.
# Now we can actually start measuring. bio3d provides the function dist.xyz() -
# but lets agree it builds character to code this ourselves.

# Consider the distance of the first (N,O) pair.
H.xyz.N[1,]
H.xyz.O[1,]

a <- H.xyz.N[1,]
b <- H.xyz.O[1,]

dist.xyz(a, b)

# or ...
sqrt( (a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2 )
# ... i.e. pythagoras theorem in 3D.


# Calculating all pairwise distances from a matrix of
# xyz coordinates sounds like a useful function.

PairDist.xyz <- function(xyzA, xyzB) {
  PD <- c()
  for (i in 1:nrow(xyzA)) {
    for (j in 1:nrow(xyzB)) {
      PD <- c(PD, dist.xyz(xyzA[i,], xyzB[j,]))
    }
  }
  return(PD)
}

# And see what we get:
D <- PairDist.xyz(H.xyz.N, H.xyz.O)
hist(D)

# let's zoom in on the shorter distances, in which we expect
# hydrogen bonds:
hist(D[D < 4.0], breaks=seq(2.0, 4.0, 0.1), xlim=c(2.0,4.0))

# There is a large peak at about 2.2Å, and another
# large peak above 3.5Å. But these are not typical hydrogen
# bond distances! Rather these are (N,O) pairs in peptide
# bonds, and within residues. That's not good, it will
# cause all sorts of problems with statistics later.
# We should exclude all distances between N of a residue
# and O of a preceeding residue, and all (N,O) pairs in the
# same residue. But for this, we need to store atom type
# and residue information with the coordinates. Our code
# will get a bit bulkier. It's often hard to find a good
# balance between terse generic code, and code that
# handles special cases.

# We could do two things:
# - add a column with residue information to the coordinates
# - add a column with atom type information
# ... but actually we also would need chain information, and
# then we really have almost everything that is contained in the record
# in the first place.

# This suggests a different strategy: let's rewrite our function
# getAtom() to return indices, not coordinates, and use the indices
# to extract coordinates, and THEN we can add all sorts of
# additional constraints.

# Here we set up the function with a default chain argument

getAtomIndex <- function(PDB, V_res, elety, chain="A") {
  # Function to access a bio3d pdb object, iterate over
  # a vector of residues and return a vector of indices
  # to matching atom elements. Nb. bio3d handles insert
  # and alt fields incorrectly: their values should not
  # be NA but " ",  i.e. a single blank. Therefore this
  # function does not test for "alt" and "insert".

  v <- c() 	#initialize as empty vector
  for (i in 1:length(V_res)) {
    res <- as.character(V_res[i])
    x <- which(PDB$atom[,"resno"] == res &
                 PDB$atom[,"chain"] == chain &
                 PDB$atom[,"elety"] == elety)
    v <- c(v, x)
  }
  return(v)
}

# test this ...
getAtomIndex(apses, H, "N")
getAtomIndex(apses, H, "O")

# That looks correct: O atoms should be stored three
# rows further down: the sequence of atoms in a PDB
# file is usually N, CA, C, O ... followed by the side
# chain coordinates.

# Now to extract the coordinates and calculate distances.
# Our function needs to take the PDB object and
# two vectors of atom indices as argument, and return
# a vector of pair-distances.

PairDist <- function(PDB, a, b) {
  PD <- c()
  for (i in 1:length(a)) {
    p <- as.numeric(PDB$atom[a[i], c("x", "y", "z")])
    for (j in 1:length(b)) {
      q <- as.numeric(PDB$atom[b[j], c("x", "y", "z")])
      PD <- c(PD, dist.xyz(p, q))
    }
  }
  return(PD)
}

# Let's see if this looks correct:

H.N <- getAtomIndex(apses, H, "N")
H.O <- getAtomIndex(apses, H, "O")
X <- PairDist(apses, H.N, H.O)
hist(X[X < 4.0], breaks=seq(2.0, 4.0, 0.1), xlim=c(2.0,4.0))

# Now we are back where we started out from, but with
# a different logic of the function that we can easily
# modify to exclude (N_i, O_i-1) distances (peptide bond),
# and (N_i, O_i) distances (within residue).

HB <- function(PDB, a, b) {
  HBcutoff <- 4.0
  PD <- c()
  for (i in 1:length(a)) {
    p <- as.numeric(PDB$atom[a[i], c("x", "y", "z")])
    res_i <- as.numeric(PDB$atom[a[i], "resno"])
    for (j in 1:length(b)) {
      q <- as.numeric(PDB$atom[b[j], c("x", "y", "z")])
      res_j <- as.numeric(PDB$atom[a[j], "resno"])
      if (res_i != res_j+1 &
          res_i != res_j     ) {
        d <- dist.xyz(p, q)
        if (d < HBcutoff) {
          PD <- c(PD, d)
        }
      }
    }
  }
  return(PD)
}

# test this:
X <- HB(apses, H.N, H.O)
hist(X)

# ... and this looks much more like the distribution we are
# seeking.

# Why did we go along this detour for coding the
# function? The point is that there are usually
# several ways to use or define datastructures and
# functions. Which one is the best way may not be
# obvious until you understand the problem better.
# At first, we wrote a very generic function that
# extracts atom coordinates to be able to compute
# with them. This is simple and elegant. But we
# recognized limitations in that we could not
# make more sophisticated selections that we needed
# to reflect our biological idea of hydrogen
# bonds. Thus we changed our datastructure
# and functions to accomodate our new requirements
# better. You have to be flexible and able to look
# at a task from different angles to succeed.

# Finally we can calculate alpha- and beta- structure
# bonds and compare them. In this section we'll explore
# different options for histogram plots.

H.N <- getAtomIndex(apses, H, "N")
H.O <- getAtomIndex(apses, H, "O")
dH <- HB(apses, H.N, H.O)

E.N <- getAtomIndex(apses, E, "N")
E.O <- getAtomIndex(apses, E, "O")
dE <- HB(apses, E.N, E.O)

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


# add color:
hist(dH, col="#DD0055")
hist(dE, col="#00AA70")



# For better comparison, plot both in the
# same window:

hist(dH, col="#DD0055")
hist(dE, col="#00AA70", add=TRUE)

# ... oops, we dind't reset the graphics parameters.
# You can either close the window, a new window
# will open with default parameters, or ...
par(opar)      # ... reset the graphics parameters

hist(dH, col="#DD0055")
hist(dE, col="#00AA70", add=TRUE)

# We see that the leftmost column of the sheet bonds
# overlaps the helix bonds. Not good. But we
# can make the colors transparent! We just need to
# add a fourth set of two hexadecimal-numbers to
# the #RRGGBB triplet. Lets use 2/3 transparent,
# in hexadecimal, 1/3 of 256 is x55 - i.e. an
# RGB triplet specied as #RRGGBB55 is only 33%
# opaque:

hist(dH, col="#DD005555")
hist(dE, col="#00AA7055", add=TRUE)

# To finalize the plots, let's do two more things:
# Explicitly define the breaks, to make sure they
# match up - otherwise they would not need to...
# see for example:

hist(dH, col="#DD005555")
hist(dE[dE < 3], col="#00AA7055", add=TRUE)

# Breaks are a parameter in hist() that can
# either be a scalar, to define how many columns
# you want, or a vector, that defines the actual
# breakpoints.
brk=seq(2.4, 4.0, 0.1)

hist(dH, col="#DD005555", breaks=brk)
hist(dE, col="#00AA7055", breaks=brk, add=TRUE)

# The last thing to do is to think about rescaling the plot.
# You notice that the y-axis is scaled in absolute frequency.
# That gives us some impression of the relative frequency,
# but it is of course skewed by observing relatively more
# or less of one type of secondary structure in a protein.
# As part of the hist() function we can rescale the values so
# that the sum over all is one: set the prameter freq=FALSE.

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
# ===========================================================
# With all the functions we have defined,
# it is easy to try this with a larger protein.
# 3ugj for example is VERY large. The calculation will take a few
# minutes:

pdb <- read.pdb("3ugj")

H <- getSecondary(pdb$helix)
E <- getSecondary(pdb$sheet)

H.N <- getAtomIndex(pdb, H, "N")
H.O <- getAtomIndex(pdb, H, "O")
dH <- HB(pdb, H.N, H.O)

E.N <- getAtomIndex(pdb, E, "N")
E.O <- getAtomIndex(pdb, E, "O")
dE <- HB(pdb, E.N, E.O)

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

# It looks more and more that the distribution is
# indeed different. Our sample is large, but derives
# from a single protein.
# To do database scale statistics, we should look
# at many more proteins. To give you a sense of how,
# let's do this for just ten proteins, taken from
# the architecture level of the CATH database for
# mixed alpha-beta proteins (see:
# http://www.cathdb.info/browse/browse_hierarchy_tree):

PDBarchitectures <- c("3A4R", "A")
names(PDBarchitectures) <- c("ID", "chain")
PDBarchitectures <- rbind(PDBarchitectures, c("1EWF","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("2VXN","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("1I3K","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("1C0P","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("3QVP","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("1J5U","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("2IMH","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("3NVS","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("1UD9","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("1XKN","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("1OZN","A"))
PDBarchitectures <- rbind(PDBarchitectures, c("2DKJ","A"))

dH <- c()
dE <- c()

for (i in 1:nrow(PDBarchitectures)) {
  pdb <- read.pdb(PDBarchitectures[i,1])
  chain <- PDBarchitectures[i,2]
  H <- getSecondary(pdb$helix)
  H.N <- getAtomIndex(pdb, H, "N", chain)
  H.O <- getAtomIndex(pdb, H, "O", chain)
  dH <- c(dH, HB(pdb, H.N, H.O))

  E <- getSecondary(pdb$sheet)
  E.N <- getAtomIndex(pdb, E, "N", chain)
  E.O <- getAtomIndex(pdb, E, "O", chain)
  dE <- c(dE, HB(pdb, E.N, E.O))
}

brk=seq(2.0, 4.0, 0.1)

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

# Why do you think these distributions are different?
# At what distance do you think H-bonds have the lowest energy?




# = 1 Tasks




# [END]
