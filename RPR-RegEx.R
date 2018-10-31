# RPR-RegEx.R
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the RPR-RegEx unit
#
# Version: 0.1
#
# Date:    2017  08  25
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 0.1    First code
#
# TODO:
#
#
# == HOW TO WORK WITH LEARNING UNIT FILES ======================================
#
# DO NOT SIMPLY  source()  THESE FILES!
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
#  going on. That's not how it works ...
#
# ==============================================================================

#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                          Line
#TOC> ----------------------------------------------
#TOC>   1        A regex example                  44
#TOC>   2        Counting lines                  111
#TOC>   2.1      Counting C-alpha atoms only     128
#TOC>   3        Code Solutions                  144
#TOC>   3.1      Counting atoms                  146
#TOC>   3.2      Counting C-alpha records        160
#TOC>
#TOC> ==========================================================================
#TOC>
#TOC>
#TOC>
#TOC>


# =    1  A regex example  =====================================================

# The canonical FASTA version of yeast Mbp1 at Uniprot
s <- ">sp|P39678|MBP1_YEAST Transcription factor MBP1 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=MBP1 PE=1 SV=1
MSNQIYSARYSGVDVYEFIHSTGSIMKRKKDDWVNATHILKAANFAKAKRTRILEKEVLK
ETHEKVQGGFGKYQGTWVPLNIAKQLAEKFSVYDQLKPLFDFTQTDGSASPPPAPKHHHA
SKVDRKKAIRSASTSAIMETKRNNKKAEENQFQSSKILGNPTAAPRKRGRPVGSTRGSRR
KLGVNLQRSQSDMGFPRPAIPNSSISTTQLPSIRSTMGPQSPTLGILEEERHDSRQQQPQ
QNNSAQFKEIDLEDGLSSDVEPSQQLQQVFNQNTGFVPQQQSSLIQTQQTESMATSVSSS
PSLPTSPGDFADSNPFEERFPGGGTSPIISMIPRYPVTSRPQTSDINDKVNKYLSKLVDY
FISNEMKSNKSLPQVLLHPPPHSAPYIDAPIDPELHTAFHWACSMGNLPIAEALYEAGTS
IRSTNSQGQTPLMRSSLFHNSYTRRTFPRIFQLLHETVFDIDSQSQTVIHHIVKRKSTTP
SAVYYLDVVLSKIKDFSPQYRIELLLNTQDKNGDTALHIASKNGDVVFFNTLVKMGALTT
ISNKEGLTANEIMNQQYEQMMIQNGTNQHVNSSNTDLNIHVNTNNIETKNDVNSMVIMSP
VSPSDYITYPSQIATNISRNIPNVVNSMKQMASIYNDLHEQHDNEIKSLQKTLKSISKTK
IQVSLKTLEVLKESSKDENGEAQTNDDFEILSRLQEQNTKKLRKRLIRYKRLIKQKLEYR
QTVLLNKLIEDETQATTNNTVEKDNNTLERLELAQELTMLQLQRKNKLSSLVKKFEDNAK
IHKYRRIIREGTEMNIEEVDSSLDVILQTLIANNNKNKGAEQIITISNANSHA"

nchar(s)
# Must be 969

# Task: Fetch the Uniprot ID by retrieving the first string that appears between
# two vertical bars ("pipes") in the header record.
#

# Develop the regular expression:
                      # Just five characters returned, so we know we are using
patt <- "^>(.{5})"    # the right functions
regmatches(s, regexec(patt, s, perl = TRUE))[[1]][2]

patt <- "^>(.*)|"    # everything to the pipe character
regmatches(s, regexec(patt, s, perl = TRUE))[[1]][2]

# Ooops - "|" is a metacharacter - we must escape it

patt <- "^>(.*)\|"    # using "\|"
# Ooops - that's not how we escape: must double the \ to send a literal
# "\" plus the character "|" to the regex engine.

patt <- "^>(.*)\\|"    # using "\\|"
regmatches(s, regexec(patt, s, perl = TRUE))[[1]][2]

# Good. Now let's first match everything that is not a "|", then match a "|"
patt <- "^>([^|]*)\\|"
regmatches(s, regexec(patt, s, perl = TRUE))[[1]][2]

# the same thing again, but capture the second match. And insist that there
# must be at least one character captured

patt <- "^>[^|]*\\|([^|]+)\\|"
# Analyze this pattern:
#    ^           anchor the match at the beginning of the line
#    >           ">" must be the first character
#    [^|]*       all-characters-except-a-vertical-bar, 0 or more times because
#                  we don't know what other versions of the string "sp"
#                  might appear. Note that within the brackets "|" is NOT a
#                  metacharacter.
#    \\|         "|" character: ouside of square brackets "|" is a metacharacter
#                  and means "OR"; we need to escape it to match a literal "|".
#    (           open parenthesis: capture what comes next ...
#       [^|]+    all-characters-except-a-vertical-bar, 1 or more times
#    )           close parenthesis: stop capturing here
#    \\|           second "|" character, escaped
regmatches(s, regexec(patt, s, perl = TRUE))[[1]][2]


# =    2  Counting lines  ======================================================

# Task: Write a function that returns the number of atoms in a PDB file. Call it
#       atomCount(). Sample data is here:
myPDB <- readLines("./data/0TST.pdb")

#       Specification:
#       Read a file from its path given as the only argument.
#       Return the number of lines in that file that begin with "ATOM  "
#       or with "HETATM".

#       Try this. Solution code is at the end of this file. Don't peek.

atomCount("./data/0TST.pdb")  # must return 6



# ==   2.1  Counting C-alpha atoms only  =======================================

# Task: write a function based on the previous one that matches only CA records,
#       i.e. it can be used to count the number of amino acids. Don't get
#       fooled by calcium atoms, or the string CA appearing elsewhere.
#       cf. https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

#       Specification:
#       Read a file from its path given as the only argument.
#       Return the number of lines in that file that have a C-alpha atom.

#       Try this. Solution code is at the end of this file. Don't peek.

CAcount("./data/0TST.pdb")  # must return 1


# =    3  Code Solutions  ======================================================

# ==   3.1  Counting atoms  ====================================================

atomCount <- function(IN) {
  # count the number of atoms in a PDB formatted file
  # Parameters:
  #     IN  chr  path of the file to read
  # Value:
  #         numeric  number of lines that match "^ATOM  " or "^HETATM"
  x <- readLines(IN)
  patt <- "(^ATOM  )|(^HETATM)"
  return(length(grep(patt, x)))
}


# ==   3.2  Counting C-alpha records  ==========================================


CAcount <- function(IN) {
  # count the number of C-alpha atoms in a PDB formatted file
  # Parameters:
  #     IN  chr  path of the file to read
  # Value:
  #         numeric  number of lines that match " CA " in position 13 - 16 of
  #                  an ATOM record.
  x <- readLines(IN)
  patt <- "^ATOM  ...... CA "
  return(length(grep(patt, x)))
}



# [END]
