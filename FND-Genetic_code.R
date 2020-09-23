# tocID <- "FND-Genetic_code.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the FND-Genetic_code unit.
#
# Version:  1.2
#
# Date:     2017  10  -  2019  01
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2020 Maintenance
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.0.1  Comment on "incomplete final line" warning in FASTA
#           1.0    First live version
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
#TOC>   Section  Title                                            Line
#TOC> ----------------------------------------------------------------
#TOC>   1        Storing the genetic code                           45
#TOC>   1.1        Genetic code in Biostrings                       63
#TOC>   2        Working with the genetic code                      94
#TOC>   2.1        Translate a sequence.                           129
#TOC>   3        An alternative representation: 3D array           212
#TOC>   3.1        Print a Genetic code table                      246
#TOC>   4        Tasks                                             272
#TOC>
#TOC> ==========================================================================


# =    1  Storing the genetic code  ============================================

# The genetic code maps trinucleotide codons to amino acids. To store it, we
# need some mechanism to associate the two representations. The most
# convenient way to do that is a "named vector" which holds the amino acid
# code and assigns the codons as names to its elements.

x <- c("M", "H", "H", "*", "*", "*")
names(x) <- c("ATG", "CAC", "CAT", "TAA", "TAG", "TGA")
x

# Then we can access the vector by the codon as name, and retrieve the
# amino acid ...

x["ATG"]
x["CAC"]
x["TAA"]

# ... or the names of elements, to retrieve the codon(s)
names(x)[x == "M"]
names(x)[x == "H"]
names(x)[x == "*"]


# ==   1.1  Genetic code in Biostrings  ========================================

# Coveniently, the standard genetic code as well as its alternatives are
# available in the Bioconductor "Biostrings" package:


if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
# Package information:
#  library(help = Biostrings)       # basic information
#  browseVignettes("Biostrings")    # available vignettes
#  data(package = "Biostrings")     # available datasets


# The standard genetic code vector
Biostrings::GENETIC_CODE

# The table of genetic codes. This information corresponds to this page
# at the NCBI:
# https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
Biostrings::GENETIC_CODE_TABLE

# Most of the alternative codes are mitochondrial codes. The id of the
# Alternative Yeast Nuclear code is "12"
Biostrings::getGeneticCode("12")  # Alternative Yeast Nuclear


# =    2  Working with the genetic code  =======================================

# We'll use Biostrings::GENETIC_CODE a lot in this script, so we'll assign it
# to a "local" variable, rather than retrieving it from the package all the
# time.

GC <- Biostrings::GENETIC_CODE

# This is a named vector of characters ...

str(GC)

# ... which also stores the alternative initiation codons TTG and CTG in
# an attribute of the vector. (Alternative initiation codons sometimes are
# used instead of ATG to intiate translation, if translation is not initiated
# at ATG thses are still translated with fMet.)

attr(GC, "alt_init_codons")

# But the key to use this vector is in the "names" which we use for subsetting
# the list of amino acids in whatever way we need.
names(GC)

# The translation of "TGG" ...
GC["TGG"]

# All stop codons
names(GC)[GC == "*"]

# All start codons
names(GC)[GC == "M"] # ... or
c(names(GC)[GC == "M"],
  attr(GC, "alt_init_codons"))


# ==   2.1  Translate a sequence.  =============================================


# I have provided a gene sequence in the data directory:
# S288C_YDL056W_MBP1_coding.fsa is the yeast Mbp1 FASTA sequence.

# read it
mbp1 <- readLines("./data/S288C_YDL056W_MBP1_coding.fsa")

# You will notice that this generates a Warning message:
#      Warning message:
#        In readLines("./data/S288C_YDL056W_MBP1_coding.fsa") :
#        incomplete final line found on './data/S288C_YDL056W_MBP1_coding.fsa'

# The reason for this is that the last character of the file is the letter "A"
# and not a "\n" line break. This file is exactly how it was sent from the
# NCBI server; I think good, defensive programming practice would have been to
# include some kind of an end-marker in the file, like a final "\n". This helps
# us recognize an incomplete transmission. Let's parse the actual sequence from
# the file, and then check for completeness.


head(mbp1)

# drop the first line (header)
mbp1 <- mbp1[-1]
head(mbp1)

# concatenate it all to a single string
mbp1 <- paste(mbp1, sep = "", collapse = "")

# how long is it?
nchar(mbp1)

# how many codons?
nchar(mbp1)/3

# That looks correct for the 833 aa sequence plus 1 stop codon. This gives us a
# first verification that the file we read is complete, the nucleotides of a
# complete ORF should be divisible by 3.

# Extract the codons. There are many ways to split a long string into chunks
# of three characters. Here we use the Biostrings  codons()  function. codons()
# requires an object of type DNAstring - a special kind of string with
# attributes that are useful for Biostrings. Thus we convert the sequence first
# with DNAstring(), then split it up, then convert it into a plain
# character vector.
mbp1Codons <- as.character(Biostrings::codons(Biostrings::DNAString(mbp1)))

head(mbp1Codons)

# now translate each codon

mbp1AA <- character(834)
for (i in seq_along(mbp1Codons)) {
  mbp1AA[i] <- GC[mbp1Codons[i]]
}

head(mbp1Codons)
head(mbp1AA)

tail(mbp1Codons)
tail(mbp1AA) # Note the stop!

# The TAA "ochre" stop codon is our second verification that the nucleotide
# sequence is complete: a stop codon can't appear internally in an ORF.

# We can work with the mbp1AA vector, for example to tabulate the
# amino acid frequencies:
table(mbp1AA)
sort(table(mbp1AA), decreasing = TRUE)

# Or we can paste all elements together into a single string. But let's remove
# the stop, it's not actually a part of the sequence. To remove the last element
# of a vector, re-assign it with a vector minus the index of the last element:
mbp1AA <- mbp1AA[-(length(mbp1AA))]
tail(mbp1AA) # Note the stop is gone!

# paste it together, collapsing the elements using an empty string as the
# separation-character (i.e.: nothing)
(Mbp1 <- paste(mbp1AA, sep = "", collapse = ""))


# =    3  An alternative representation: 3D array  =============================


# We don't use 3D arrays often - usually just 2D tables and data frames, so
# here is a good opportunity to review the syntax of 3D arrays with a
# genetic code cube:

# Initialize, using A G C T as the names of the elements in each dimension
cCube <- array(data     = character(64),
               dim      = c(4, 4, 4),
               dimnames = list(c("A", "G", "C", "T"),
                               c("A", "G", "C", "T"),
                               c("A", "G", "C", "T")))

# fill it with amino acid codes using three nested loops
for (i in 1:4) {
  for (j in 1:4) {
    for (k in 1:4) {
      myCodon <- paste(dimnames(cCube)[[1]][i],
                       dimnames(cCube)[[2]][j],
                       dimnames(cCube)[[3]][k],
                       sep = "",
                       collapse = "")
      cCube[i, j, k] <- GC[myCodon]
    }
  }
}

# confirm
cCube["A", "T", "G"] # methionine
cCube["T", "T", "T"] # phenylalanine
cCube["T", "A", "G"] # stop (amber)



# ==   3.1  Print a Genetic code table  ========================================


# The data structure of our cCube is well suited to print a table. In the
# "standard" way to print the genetic code, we write codons with the same
# second nucleotide in columns, and arrange rows in blocks of same
# first nucleotide, varying the third nucleotide fastest. This maximizes the
# similarity of adjacent amino acids in the table if we print the
# nucleotides in the order T C A G. It's immidiately obvious that the code
# is not random: the universal genetic code is exceptionally error tolerant in
# the sense that mutations (or single-nucleotide translation errors) are likely
# to result in an amino acid with similar biophysical properties as the
# original.

nuc <- c("T", "C", "A", "G")

# (calling variables f, s, t to indicate first, second, and third position ...)
for (f in nuc) {      # first varies in blocks
  for (t in nuc) {    # third varies in columns
    for (s in nuc) {  # second varies in rows
      cat(sprintf("%s%s%s: %s   ", f, s, t, cCube[f, s, t]))
    }
    cat("\n")
  }
  cat("\n")
}


# =    4  Tasks  ===============================================================


# Task: What do you need to change to print the table with U instead
#         of T? Try it.


# Task: Point mutations are more often transitions (purine -> purine;
#         pyrimidine -> pyrimidine) than transversions (purine -> pyrimidine;
#         pyrimidine -> purine), even though twice as many transversions
#         are possible in the code. This is most likely due a deamination /
#         tautomerization process that favours C -> T changes. If the code
#         indeed minimizes the effect of mutations, you would expect that
#         codons that differ by a transition code for more similar amino acids
#         than codons that differ by a transversion. Is that true? List the set
#         of all amino acid pairs that are encoded by codons with a C -> T
#         transition. Then list the set of amino acid pairs with a C -> A
#         transversion. Which set of pairs is more similar?


# Task: How many stop codons do the two mbp1-gene derived amino acid sequences
#         have if you translate them in the 2. or the 3. frame?


# Task: How does the amino acid composition change if you translate the mbp1
#         gene with the Alternative Yeast Nuclear code that is used by the
#         "GTC clade" of fungi?
#         (cf. https://en.wikipedia.org/wiki/Alternative_yeast_nuclear_code )

# Solution:

    # Fetch the code
    Biostrings::GENETIC_CODE_TABLE
    Biostrings::GENETIC_CODE_TABLE$name[Biostrings::GENETIC_CODE_TABLE$id=="12"]
    altYcode <- Biostrings::getGeneticCode("12")

    # what's the difference?
    (delta <- which(Biostrings::GENETIC_CODE != altYcode))

    Biostrings::GENETIC_CODE[delta]
    altYcode[delta]

    # translate
    altYAA <- character(834)
    for (i in seq_along(mbp1Codons)) {
      altYAA[i] <- altYcode[mbp1Codons[i]]
    }

    table(mbp1AA)
    table(altYAA)

# Task: The genetic code has significant redundacy, i.e. there are up to six
#         codons that code for the same amino acid. Write code that lists how
#         many amino acids are present how often i.e. it should tell you that
#         two amino acids are encoded only with a single codon, three amino
#         acids have six codons, etc. Solution below, but don't peek. There
#         are many possible ways to do this.
#
#
# Solution:
( x <- table(table(Biostrings::GENETIC_CODE)) )

# confirm
sum(x * as.numeric(names(x)))



# [END]
