# RPR-FASTA.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-FASTA unit.
#
# Version:  1.0
#
# Date:     2017  10  14
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    New unit.
#
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
#TOC>   Section  Title                 Line
#TOC> -------------------------------------
#TOC>   1        Reading FASTA           39
#TOC>   2        Interpreting FASTA     227
#TOC>   3        Writing FASTA          248
#TOC> 
#TOC> ==========================================================================
 



# =    1  Reading FASTA  =======================================================

# FASTA is a text based format, structured in lines that are separated by
# line-feed or paragraph-break characters. Which one of these is used, depends
# on your operating system. But Rs readLines() function knows how to handle
# these correctly, accross platforms. Don't try to read such files "by hand".
# Here is the yeast Mbp1 gene, via SGD.

file.show("./data/S288C_YDL056W_MBP1_coding.fsa")
myFASTA <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")

# The warning is generated because the programmer who implemented the code to
# write this FASTA file neglected to place a line-break character after the last
# sequence character.

head(myFASTA)

# Note that there are NO line-break characters ("\n") at the end of these
# strings, readLines() has "consumed" them while reading.

tail(myFASTA)

# Also note that the last line has fewer characters - this means readLines()
# imported the whole line, despite it not being terminated.

# It's very straightforward to work with such data, for example by collapsing
# everything after the first line into a single string ...

f <- c(myFASTA[1], paste(myFASTA[-1], sep = "", collapse = ""))

f[1]
nchar(f[2])

# ... but this is making assumptions that everything in line 2 until the end IS
# sequence, the whole sequence and nothing but sequence. That assumption can
# break down in many ways:
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
# assumptions that you are making. Here is the structure of a FASTA file,
# specified with as few assumptions as possible.
#
#  (1) it contains characters;
#  (2) there might be lines that begin with characters other than
#         ">", these should be discarded;
#  (3) it contains one or more consecutive lines that are sequence blocks;
#  (4) each sequence block has one or more header lines;
#  (5) header lines start with ">";
#  (6) no actual sequence data begins with a ">";
#  (7) header lines can contain any character;
#  (8) sequence lines only contain letters, "-" (gap characters), or "*" (stop).
#
# This suggests to parse as follows:
# - drop all lines that don't begin with ">" or a letter
# - identify consecutive lines that begin ">" and consecutive lines
#     that do not begin ">"
# - collapse each set of consecutive lines in-place
# - drop all remaining lines. In this result the odd-indexed elements
#     are headers, and the even-indexed elements are sequences.

# Let's code this as a function. We need some tool that identifies consecutive
# lines of something. The rle() (run-length encoding) function does this. It
# returns a vector of the length of "runs" in its input:

myPets <- c("ant", "bat", "bat", "bat", "cat", "ant", "ant")
(runs <- rle(myPets))

# The cumsum() (cumulative sum) function turns these numbers into indices
# on our original vector.

(idx <- cumsum(runs$lengths))
myPets[idx]   # note that this is NOT unique ... "ant" appears twice, because
              # there were two separate runs of ants in our input.

# So far so good. But our FASTA file's lines are ALL different, so all the runs
# will only have length 1 ...

rle(myFASTA)$lengths

# How do we deal with that? Obviously we need to actually analyze the strings we
# are working with. grepl(<pattern>, <x>) is exactly what we need here. It
# produces a vector of booleans, of the same length as the input vector <x>,
# which is TRUE if the element matches the <pattern>, FALSE if not.

grepl("^>", myFASTA)  # "^>" is a regular expression that means: ">" at the
                      # beginning ("^") of the line.

(runs <- rle(grepl("^>", myFASTA)))

# Translating that into start positions of blocks takes a bit of bookkeeping:
# the first start has index 1, the following starts can be calculated from
# cumsum()'s and $length's.
(starts <- c(1, (cumsum(runs$lengths)[-length(runs$lengths)] + 1)))

# ... and with that, we can parse our FASTA data. We take the specification
# above and translate it into code. That's how we develop code: write up step by
# instructions as comments, then implement them one by one.

# Here is an example
FA <- c(">head1 part a", ">head1 part b", "abcdef", "ghi", # two headers
        "",                                                # empty line
        ">head2", "jkl",                                   # one header
        ">head3", "mno", "pqrs")                           # two sequence lines

# - drop all lines that don't begin with ">" or a letter, "-", or "*"
FA <- FA[grepl("^[A-Za-z>*-]", FA)]

# - identify consecutive lines that begin ">" and consecutive lines
#     that do not begin ">"
runs <- rle(grepl("^>", FA))
starts <- c(1, (cumsum(runs$lengths)[-length(runs$lengths)] + 1))

# - collapse each set of consecutive lines in-place

for (i in seq_along(starts)) {
  FA[starts[i]] <- paste(FA[starts[i]:(starts[i] + runs$lengths[i] - 1)],
                         sep ="",
                         collapse = "")
}

# - drop all remaining lines.
FA <- FA[starts]

# In this resulting vector the odd-indexed elements
#     are headers, and the even-indexed elements are sequences.

# As a function:

readFASTA <- function(IN) {
  # Read a FASTA formatted file from IN, remove all non-header, non-sequence
  # element, return collapsed sequences.
  # Parameters:
  #    IN  chr   Input file name (or connection)
  # Value:
  #    chr vector  in which the odd-indexed elements are headers, and the
  #                even-indexed elements are sequences.

  FA <- readLines(IN)
  FA <- FA[grepl("^[A-Za-z>*-]", FA)]

  runs <- rle(grepl("^>", FA))
  starts <- c(1, (cumsum(runs$lengths)[-length(runs$lengths)] + 1))

  for (i in seq_along(starts)) { # collapse runs in-place
    FA[starts[i]] <- paste(FA[starts[i]:(starts[i] + runs$lengths[i] - 1)],
                           sep ="",
                           collapse = "")
  }

  # return collapsed lines
 return(FA[starts])
}

# Try this: Let's try to use only the first 3 elements of myFASTA ... it's a
# lengthy sequence. But how? We don't have a file with that contents and the
# function expects to read from a file. Do we need to write myFASTA[1:3] to a
# temporary file and then read it? We could - but wherever a file is expected we
# can also pass in a "text connection" from an object in memory, with the
# textConnection() function, like so:

readFASTA(textConnection(myFASTA[1:3]))

# Here is a "real" example - a multi FASTA file of aligned APSES domain
# sequences:

(refAPSES <- readFASTA("./data/refAPSES.mfa"))

# Subset all headers:
refAPSES[seq(1, length(refAPSES), by = 2)]

# Subset the sequence with "P39678" in the header
refAPSES[grep("P39678", refAPSES) + 1]  # grep() the string and add 1



# =    2  Interpreting FASTA  ==================================================


# FASTA files are straightforward to interpret - just one thing may be of note:
# when working with strings, we can use substr(<string>, <start>, <stop>) to
# extract substrings, but more often we expand the string into a vector of
# single characters with strsplit(<string>, ""). strsplit() returns a list,
# to accommodate that <string> could be a vector of many elements, therafore
# we usually unlist() the result if we use it only on a single string.

# Example: How many positive charged residues in "MBP1_SACCE"?

s <- unlist(strsplit(refAPSES[grep("MBP1_SACCE", refAPSES) + 1], ""))
head(s)
sum(grepl("[HKR]", s)) # 20 (+) charged residues. grepl() returns TRUE and FALSE
                       # for the characters, sum() coerces to 1 and 0
                       # respectively, and that gives us the result.

100 * sum(grepl("[HKR]", s)) / length(s) # in percent: 20.2 %


# =    3  Writing FASTA  =======================================================


# Writing FASTA files mostly just the revrese reverse of reading, with one
# twist: we need to break the long sequence string into chunks of the desired
# width. The FASTA specification calls for a maximum of 120 characters per line,
# but writing out much less than that is common since it allows to comfortably
# view lines on the console, or printing them on a sheet of paper (do we still
# do that actually?). How do we break a string into chunks? A combination of
# seq(<from>, <to>, <by>) with substring(<string>, <start>, <stop>) will work
# nicely. (Note that substring() is vectorized, whereas substr() is not!) As we
# loop through our FASTA object in memory, we can build the output by c()'ing
# blocks of header + sequence to each other. For VERY large objects this might
# be slow - in that case, we might want to precalculate the size of the output
# object. But that's more of a hypothetical consideration.

s <- refAPSES[2]
nchar(s)
w <- 30     # width of chunk
(starts <- seq(1, nchar(s), by = w))      # starting index of chunk
(ends <- c((starts - 1)[-1], nchar(s)))   # ending index of chunk

# Task: Is this safe? What happens if nchar(s) is shorter than w?
#       What happens if nchar(s) is an exact multiple of w?

substring(s, starts, ends)

# Here's the function ...

writeFASTA <- function(s, OUT = stdout(), width = 60) {
  # Write an object "s" that contains one or more header/sequence pairs to file.
  # Parameters:
  #    s      chr   Vector with a FASTA header string in odd elements,
  #                    sequence in one-letter code in even elements.
  #    OUT    chr   connection to be written to; defaults to stdout() i.e.
  #                 output is written console.
  #    width  int   max number of sequence characters per line of output.
  # Value:
  #           NA    Invoked for side effect of writing data to file

  txt <- character()
  idx <- seq(1, length(s), by = 2)
  for (i in idx) {
    txt <- c(txt, s[i])                              # add header line to txt
    starts <- seq(1, nchar(s[i + 1]), by = width)    # starting indices of chunks
    ends <- c((starts - 1)[-1], nchar(s[i + 1]))     # ending indices of chunks
    txt <- c(txt, substring(s[i + 1], starts, ends)) # add chunks to txt
  }
  writeLines(txt, OUT)

}

# Let's try this. We don't define OUT, so the result is written to the console
# by default. Defualt width for sequence is 60 characters

writeFASTA(refAPSES)



# [END]
