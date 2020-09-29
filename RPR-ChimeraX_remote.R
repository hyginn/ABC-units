# tocID <- "RPR-ChimeraX_remote.R"
#
# Purpose:  A Bioinformatics Course:
#              R code demonstrating remote scripting of ChimeraX.
#
# Version:  1.0
#
# Date:     2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    First ABC units version
#
#
# TODO:
#    %-encode and escape quotes, or just pass-through?
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
#TOC>   Section  Title                                  Line
#TOC> ------------------------------------------------------
#TOC>   1        ChimeraX REMOTE SCRIPTING                40
#TOC>   1.1        Defining a Port                        58
#TOC>   1.2        Open ChimeraX                          80
#TOC>   2        WORKED EXAMPLE: SUPERPOSITION           112
#TOC> 
#TOC> ==========================================================================


# =    1  ChimeraX REMOTE SCRIPTING  ===========================================


# One of the cool features of ChimeraX is that it can be driven by Python code,
# both within a running session and through Python scripts. What I find even
# cooler though is that ChimeraX can be driven from any programming language via
# its remote control function that can listen to commands sent from any other
# application. The interface that is used here is the standard REST (method) -
# the GET and POST verbs that ubiquitously underly the communication of clients
# and servers on the Web.

# In order to establish the communication between this script and ChimeraX, all
# we need to do is:
#  - open ChimeraX;
#  - tell it to listen on a specific "port";
#  - send commands to that port via httr::


# ==   1.1  Defining a Port  ===================================================

# The httr:: package needs to be available

if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes
#  data(package = "httr")     # available datasets

# We need to think od a port. Any available port number between 49152-65535 is
# fine. We'll choose 61803 because that's the fractional part of the golden
# ratio. But one could choose another.

CXPORT <- 61803

# Check that our current version of R supports sockets (default since V 3.3)
capabilities("sockets")   # MUST be TRUE. If not, don't continue.


# ==   1.2  Open ChimeraX  =====================================================

#  - Open a fresh, new session of ChimeraX
#  - type:
#
#       remotecontrol rest start port 61803
#
#    ... or whatever the value of CXPORT is.

# Now watch what happens in ChimeraX when you execute the following line:
( x <- httr::GET("http://127.0.0.1:61803/run?command=open+1BM8") )

# The .utilities.R script includes the function CX(), based on this principle,
# through which you can send commands to ChimeraX

CX("camera sbs")
CX("lighting soft")
CX("color sequential #1 & protein target abc palette powderblue:orchid:white")

# The command echos Chimera's response if the parameter "quietly" is
# FALSE (default), and we can silence output with quietly = TRUE :
CX("info models #1 attribute num_residues")
CX("info models #1 attribute num_residues", quietly = TRUE)

# Either way, the command also returns Chimera's responses "invisibly";
# i.e. we can use the results by assigning the output to a variable:
hBonds <- CX("hbonds #1 & protein makePseudobonds false log true", quietly=TRUE)
x <- read.table(file = textConnection(hBonds), skip = 9,
                blank.lines.skip = TRUE, fill = TRUE)
hist(x[,13], main="H-bonds", xlab="D···A (Å)", ylab="counts", col="#c9dcff")


# =    2  WORKED EXAMPLE: SUPERPOSITION  =======================================

# We superimpose the 1BM8 structure with the 1DUX crystal structure to be able
# to explore possible DNA binding regions in 1BM8

# The model for 1BM8 is already open as model 1  (#1)
CX("hide #1 cartoons")        # hide model 1 cartoon representation
CX("open 1DUX")               # assume this is opened as model #2
CX("hide #2")                 # hide everything ...
CX("select #2/C")             # chain c (protein)
CX("show sel cartoons")       # ... and show cartoons of chain c (protein)
CX("color sequential sel target c palette steelblue:darkmagenta")
CX("view #2/C")               # re-center the display
CX("cofr #2/C:62@CA")         # set pivot to an interface residue
CX("select #2/A,B & nucleic-acid") # chains A, B are the cognate DNA
CX("style sel stick")
CX("show sel target ab")      # show atoms/bonds
CX("color sequential #2/A & nucleic-acid target ab palette teal:lightcyan")
CX("color sequential #2/B & nucleic-acid target ab palette teal:lightcyan")
CX("surface sel enclose sel") # compute joint accessible surface of both chains
CX("transparency 50")
CX("select clear")

# Now superimpose the 1BM8 chain onto 1DUX chain C
CX("show #1 cartoons")
CX("matchmaker #1/A to #2/C pairing ss")  # the actual superposition

# study the general layout, and the position of the 1mb8 secondary structure
# elements relative to 1DUX

# Let's examine side chain orientations in more detail
CX("hide #2/C cartoons")  # hide the 1DUX protein

# select all residues in 1BM8 that are within 3.5 A of the DNA chains (a, b)
CX("select zone #2/A,B 3.5 #1 & protein residues true")
CX("~select sel & H")  # de-select H atoms
CX("show sel target ab")
CX("size stickRadius 0.4")
CX("select clear")

# The overall architecture of the Mbp1 APSES domain is a good match for the Elk
# transcription factor binding mode; the detailed conformations of side chains
# would need to change only to a minor degree. There is a very significant
# degree of structural similarity; remarkable, given that the DNA is not the
# target sequence of the Mbp1 transcription factor, AND the 1MB8 structure was
# determined without a DNA ligand.

CX("remotecontrol rest stop")  # release the socket
# Done.



# [END]
