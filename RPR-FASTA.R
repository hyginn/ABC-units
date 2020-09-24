# tocID <- "RPR-FASTA.R"
#
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-FASTA unit.
#
# Version:  1.1
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.1    2020 Maintenance. Rewrite validation logic. Add data
#                  to utilities. Define AACOLS
#           1.0    New unit.
#
#
# TODO: Make a simple solution first, then extend it to error checking, and
#       to handle .mfa files.
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
#TOC>   Section  Title                                 Line
#TOC> -----------------------------------------------------
#TOC>   1        Reading and validating FASTA            45
#TOC>   1.1        Validating FASTA                      81
#TOC>   2        Parsing FASTA                          225
#TOC>   3        Interpreting FASTA                     245
#TOC>   4        Writing FASTA                          272
#TOC> 
#TOC> ==========================================================================


# =    1  Reading and validating FASTA  ========================================

# FASTA is a text based format, structured in lines that are separated by
# line-feed or paragraph-break characters. Which one of these is used, depends
# on your operating system. But R's readLines() function knows how to handle
# these correctly, accross platforms. Don't try to read such files "by hand".
# Here is the yeast Mbp1 gene, via SGD.

file.show("./data/S288C_YDL056W_MBP1_coding.fsa")
faMBP1 <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")

# The warning is generated because the programmer at the NCBI who implemented
# the code to write this FASTA file neglected to place a line-break character
# after the last sequence character. While this is not technically incorrect,
# it is poor practice.

head(faMBP1)

# Note that there are NO line-break characters ("\n") at the end of these
# strings, even though they were present in the original file. readLines()
# has "consumed" these characters while reading - but every single line is in
# a vector of its own.

tail(faMBP1)

# Also note that the last line has fewer characters - this means readLines()
# imported the whole line, despite it not being terminated by "\n".

# It's very straightforward to work with such data, for example by collapsing
# everything except the first line into a single string ...

f <- c(faMBP1[1], paste(faMBP1[-1], sep = "", collapse = ""))

f[1]
nchar(f[2])

# ==   1.1  Validating FASTA  ==================================================

# The code above is making the assumption that everything from line 2 until
#  the end IS sequence, the whole sequence and nothing but sequence.
#  That assumption can break down in many ways:
#
#  - there could be more than one header line. The specification says otherwise,
#       but some older files use multiple, consecutive header lines. You don't
#       want that to end up in your sequence.
#  - this could be not a FASTA file at all. It could be raw sequence, a
#       different sequence file format, or a wholly different file altogether.
#       If you look at the file, you can immediately tell, but if you are
#       reading the file in a complex workflow, your could easily import wrong
#       data into your analysis.
#  - there could be more than one sequence in the file. Such Multi-FASTA files
#       occur commonly, as downloads of ORFs from genome regions or other
#       sets of genes or proteins, or as the input / output for multiple
#       sequence alignment programs.
#
# Data "from the wild" can (and usually does) have the most unexpected
# variations and it is really, really important to be clear about the
# assumptions that you are making. It is possible to "fix" things, according
# to the "Robustness Principle" :
#      "Be conservative in what you send,
#       be liberal in what you accept".
#       (cf. https://en.wikipedia.org/wiki/Robustness_principle )
# ... but if you think about this, that's actually a really poor idea,
# which is much more likely to dilute standards, make unwarranted
# assumptions, and allow errors to pass silently and corrupt data.
#
# Let's discard this principle on the trash-heap of
# things-that-sound-like-a-good-idea-but-aren't. What we do instead is test,
# identify problems, and follow the principle: "crash early, crash often". Of
# course I can write code that would reformat any possible input as a FASTA
# file - but what good will it do me if it parses the file I receive
# from a server into FASTA format like:
#
#   >404- Page Not Found</title</head>
#   dyh-PagentfndhpThepageyreqesteddesnteistnthisserverCheckthe
#   spellingrcntacttheadministratrsdyhtml
#
# Therefore, we write ourselves a FASTA checker that will enforce the following:
#   (1) a FASTA file contains one or more sequences separated by zero or
#       more empty lines
#   (2) a sequence contains one header line followed by
#       one or more sequence lines
#   (3) a sequence line contains one or more uppercase or lowercase single
#       letter amino acid codes, hyphens (gap character), or * (stop).
#
#   Anything else should generate an error.

#   (Case 1): Header(s) exist
fX <- c("ABC",
        "defghi",
        "klmnpq")
sel <- grepl("^>", fX)  # "^>" is a regular expression that
                        # means: the exact character ">" at the
                        # beginning ("^") of the line.
if ( ! any(sel) ) { stop("no header lines in input.") }


#   (Case 2) No adjacent header lines
fX <- c(">ABC",
        ">123",
        "defghi",
        "klmnpq")
sel <- grepl("^>", fX)
sel <- sel[- length(sel)] & sel[-1] # comparing shifted vectors
if ( any(sel)) { stop("adjacent header lines in input.") }

#   (Case 3.1) all sequence lines contain only valid characters
#              (constants for valid characters AAVALID, NUCVALID, and NUCAMBIG
#               are defined with the .utilities.R script)
AAVALID
fX <- c(">ABC",
        "def ;-) ghi",
        "klmnpq")
myRegex <- sprintf("[^%s]", AAVALID)  # NOT a valid character
sel <- ! grepl("^>", fX)              # NOT headers
if (any(grepl(myRegex, fX[sel]))) {
  stop("invalid chracter(s) outside of header lines.")
}

#   (Case 3.2) all headers are followed directly by
#              at least one letter of sequence
fX <- c(">ABC",
        "",
        ">123",
        "defghi",
        "klmnpq")
sel <- grep("^>", fX) + 1             # indexes of headers + 1
myRegex <- sprintf("[%s]+", AAVALID)  # at least one valid character
if (! all(grepl(myRegex, fX[sel]))) {
  stop("a header has no adjacent sequence.")
}
# Ah, you might ask - couldn't we just have dropped all empty lines, and
# then caught this in Case 2? No - for two reasons: we would still miss headers
# at the end of file, and, we would have changed the line numbering - and
# ideally our "production" function will create information about where the
# error is to be found.


# Now combine this into a function ...

val <- function(fa) {

  if ( ! any(grepl("^>", fa)) ) {
    stop("no header lines in input.")
  }

  sel <- grepl("^>", fa)
  if ( any(sel[- length(sel)] & sel[-1])) {
    stop("adjacent header lines in input.")
  }

  sel <- ! grepl("^>", fa)
  if ( any(grepl(sprintf("[^%s]", AAVALID), fa[sel]))) {
    stop("invalid chracter(s) outside of header lines.")
  }

  sel <- grep("^>", fa) + 1
  if (! all(grepl(sprintf("[%s]+", AAVALID), fa[sel]))) {
    stop("a header has no adjacent sequence.")
  }

  return(invisible(NULL))
}

# Here is an example
FA <- c(">head1",
        "acdef",
        "ghi",
        "",
        ">head2",
        "kl",
        ">head3",
        "mn",
        "pqrs")
validate(FA)     # ... should not create an error


# a somewhat more elaborate validateFA() function was loaded with the
# ./utilities.R script

# =    2  Parsing FASTA  =======================================================

# Once we have validated our assumptions about our input, it's quite
# painless to parse it. I have put this together as a function and the function
# gets loaded from ./.utilities.R
#

# Lets try this:
#   - the first 3 elements of faMBP1:
readFASTA(faMBP1[1:3])

#   - a multi FASTA file of aligned APSES domain sequences:

refAPSES <- readFASTA("./data/refAPSES.mfa")

# Subset the sequence with "P39678" in the header
refAPSES[grep("P39678", refAPSES$head) ,]



# =    3  Interpreting FASTA  ==================================================


# FASTA files are straightforward to interpret - just one thing may be of note:
# when working with strings, we can use substr(<string>, <start>, <stop>) to
# extract substrings, but more often we expand the string into a vector of
# single characters with strsplit(<string>, ""). strsplit() returns a list,
# to accommodate that <string> could be a vector of many elements, therefore
# we usually unlist() the result if we use it only on a single string.

# Example: How many positive charged residues in "MBP1_SACCE"?

s <- unlist(strsplit(refAPSES$seq[grep("MBP1_SACCE", refAPSES$head)], ""))
s

sum(grepl("[HKR]", s)) # 20 (+) charged residues. grepl() returns TRUE and FALSE
                       # for the characters, sum() coerces to 1 and 0
                       # respectively, and that gives us the result.

100 * sum(grepl("[HKR]", s)) / length(s) # in percent: 20.2 %

# residue distribution
x <- factor(s, levels = names(AACOLS))
pie(table(x)[names(AACOLS)], col = AACOLS)



# =    4  Writing FASTA  =======================================================


# Writing FASTA files mostly just the revrese reverse of reading, with one
# twist: we need to break the long sequence string into chunks of the desired
# width. The FASTA specification calls for a maximum of 120 characters per line,
# but writing out much less than that is common, since it allows to comfortably
# view lines on the console, or printing them on a sheet of paper (do we still
# do that actually?). How do we break a string into chunks? A combination of
# seq(<from>, <to>, <by>) with substring(<string>, <start>, <stop>) will work
# nicely. (Note that substring() is vectorized, whereas substr() is not!) As we
# loop through our FASTA object in memory, we can build the output by c()'ing
# blocks of header + sequence to each other. For VERY large objects this might
# be slow - in that case, we might want to precalculate the size of the output
# object. But that's more of a hypothetical consideration.

( s <- refAPSES$seq[2] )
nchar(s)
w <- 30     # width of chunk
(starts <- seq(1, nchar(s), by = w))      # starting index of chunk
(ends <- c((starts - 1)[-1], nchar(s)))   # ending index of chunk

# Task: Is this safe? What happens if nchar(s) is shorter than w?
#       What happens if nchar(s) is an exact multiple of w?

substring(s, starts, ends)
# confirm that the output contains the first and last residue, and both
# residues adjacent to the breaks

# As always, the function has been defined in ".utilities.R" for to use
# any time...  type   writeFASTA  to examine it.

# Let's try this...

writeFASTA(refAPSES, width = 40)

# roundtrip for validation: write refAPSES with a different format,
# read it back in - the new dataframe must be identical
# to the original dataframe.
fname <- tempfile()
writeFASTA(refAPSES, fn = fname, width = 30)
identical(refAPSES, readFASTA(fname))

# ...works for me  :-)


# [END]
