# ABC_makeMYSPElist.R
#
# Purpose:  Create a list of genome sequenced fungi with protein annotations and
#               Mbp1 homologues.
#
# Version: 1.1.2
#
# Date:    2016 09 - 2017 09
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.1.2  Moved BLAST.R to ./scripts directory
# V 1.1    Update 2017
# V 1.0    First code 2016
#
# TODO:
#
#   type out workflow
#
# ==============================================================================
#
# DO NOT  source()  THIS FILE!
#
# This file is code I provide for your deeper understanding of a process and
# to provide you with useful sample code. It is not actually necessary for
# you to run this code, but I encourage you to read it carefully and discuss
# if there are parts you don't understand.
#
# Run the commands that interact with the NCBI servers only if you want to
# experiment specifically with the code and/or parameters. I have commented out
# those parts. If you only want to study the general workflow, just load()
# the respective intermediate results.
#

#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                               Line
#TOC> ---------------------------------------------------
#TOC>   1        The strategy                          54
#TOC>   2        GOLD species                          66
#TOC>   2.1      Initialize                            71
#TOC>   2.2      Import                                77
#TOC>   2.3      Unique species                       129
#TOC>   3        BLAST species                        171
#TOC>   3.1      find homologous proteins             178
#TOC>   3.2      Identify species in "hits"           202
#TOC>   4        Intersect GOLD and BLAST species     247
#TOC>   5        Cleanup and finish                   265
#TOC>
#TOC> ==========================================================================


#TOC>
#TOC>

# =    1  The strategy  ========================================================

# This script will create a list of "MYSPE" species and save it in an R object
# MYSPEspecies that is stored in the data subdirectory of this project from where
# it can be loaded. The strategy is as follows: we download a list of all
# genome projects and then select species for which protein annotations are
# available - i.e. these are all genome-sequenced species that have been
# annotated. Then we search for fungal species that have homologues to MBP1.
# Then we intersect the two lists to give us genome-sequenced species that
# also have Mbp1 homologues ...


# =    2  GOLD species  ========================================================

#  Fetch and parse the Genomes OnLine Database of the Joint Genome Institute
#  (https://gold.jgi.doe.gov/). Use the data that is hosted at the NCBI.

# ==   2.1  Initialize  ========================================================
if (! require(httr)) { # httr provides interfaces to Webservers on the Internet
    install.packages("httr")
    library(httr)
}

# ==   2.2  Import  ============================================================

# The URL of the genome data directory at the NCBI:
# is https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS
# Note the relative size of the prokaryotes and the eukaryotes data.

# What's in this directory?
URL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/README"
GOLDreadme <- readLines(URL) # read the file into a vector
cat(GOLDreadme, sep = "\n")  # display the contents

# Retrieve the file "eukaryotes" via ftp from the NCBI ftp server and put it
# into a dataframe. This will take a few moments.
# URL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"
# GOLDdata <- read.csv(URL,
#                      header = TRUE,
#                      sep = "\t",
#                      stringsAsFactors = FALSE)
# save(GOLDdata, file="data/GOLDdata.RData")
# or ...
load(file="data/GOLDdata.RData")


# What columns does the table have, how is it structured?
str(GOLDdata)

# What groups of organisms are in the table? How many of each?
table(GOLDdata$Group)

# What subgroups of fungi do we have?
table(GOLDdata$SubGroup[GOLDdata$Group == "Fungi"])

# How many of the fungi have protein annotations? The README file told us that
# the column "Proteins" contains "Number of Proteins annotated in the assembly".
# Looking at a few ...
head(GOLDdata$Proteins, 30)
# ... we see that the number varies, and some have a hyphen, i.e. no
# annotations. The hyphens make this a char type column (as per: all elements
# of a vector must have the same type). Therefore we can't read this as numbers
# and filter by some value > 0. But we can filter for all genomes that don't
# have the hyphen:
sum(GOLDdata$Proteins[GOLDdata$Group == "Fungi"] != "-")

# Subset the data, with fungi that have protein annotations
GOLDfungi <- GOLDdata[GOLDdata$Group == "Fungi" &
                          GOLDdata$Proteins != "-" , ]

# check what we have in the table
nrow(GOLDfungi)
head(GOLDfungi)


# ==   2.3  Unique species  ====================================================

# For our purpose of defining species, we will select only species, not strains
# from this list. To do this, we pick the first two words i.e. the systematic
# binomial name from the "X.Organism.Name" column, and then we remove redundant
# species. Here is a function:
#

getBinom <- function(s) {
    # Fetch the first two words from a string.
    # Parameters:
    #   s: char  a string which is expected to contain a binomial species name
    #            as the first two words, possibly followed by other text.
    # Value: char  the first two words separated by a single blank
    #
    x <- unlist(strsplit(s, "\\s+"))     # split s on one or more whitespace
    return(paste(x[1:2], collapse=" "))  # return first two elements
}

# iterate through GOLDdata and extract species names
GOLDspecies <- character()
for (i in 1:nrow(GOLDfungi)) {
    GOLDspecies[i] <- getBinom(GOLDfungi$X.Organism.Name[i])
}
head(GOLDspecies)
length(GOLDspecies)

# N.b. this would be more efficiently (but perhaps less explicitly) coded with
# one of the apply() functions, instead of a for-loop.
# GOLDspecies <- unlist(lapply(GOLDfungi$X.Organism.Name, getBinom))

# Species of great interest may appear more than once, one for each sequenced
# strain: e.g. brewer's yeast:
sum(GOLDspecies == "Saccharomyces cerevisiae")

# Therefore we use the function unique() to throw out duplicates. Simple:
GOLDspecies <- unique(GOLDspecies)

length(GOLDspecies)
# i.e. we got rid of about 40% of the species by removing duplicates.


# =    3  BLAST species  =======================================================
#
# Next, we filter our list by species that have homologues to the yeast Mbp1
# gene. To do this we run a BLAST search to find all related proteins in any
# fungus. We list the species that appear in that list, and then we select those
# that appear in our GOLD table as well.
#
# ==   3.1  find homologous proteins  ==========================================
#
# Use BLAST to fetch proteins related to Mbp1 and identify the species that
# contain them.

# Scripting against NCBI APIs is not exactly enjoyable - there is usually a fair
# amount of error handling involved that is not supported by the API in a
# principled way but requires rather ad hoc solutions. The code I threw together
# to make a BLAST interface (demo-quality, not research-quality) is in the file
# ./scripts/BLAST.R Feel encouraged to study how this works. It's a pretty
# standard task of communicating with servers and parsing responses - everyday
# fare in thebioinformatics lab. Surprisingly, there seems to be no good BLAST
# parser in currently available packages.

# source("./scripts/BLAST.R")   # load the function and its utilities
# Use BLAST() to find yeast Mbp1 homologues in other fungi in refseq
# BLASThits <- BLAST("NP_010227",                  # Yeast Mbp1 RefSeq ID
#                    db = "refseq_protein",        # database to search in
#                    nHits = 3000,                 # 720 hits in 2017
#                    E = 0.01,                     #
#                    limits = "txid4751[ORGN]")    # = fungi
# save(BLASThits, file="data/BLASThits.RData")
load(file="data/BLASThits.RData")

# ==   3.2  Identify species in "hits"  ========================================

# This is a very big list that can't be usefully analyzed manually. Here
# we are only interested in the species names that it contains.

# How many hits in the list?
length(BLASThits$hits)

# Let's look at a hit somewhere down the list
str(BLASThits$hit[[277]])

# A fair amount of parsing has gone into the BLAST.R code to prepare the results
# in a useful way. The species information is in the $species element of every
# hit.

# Run a loop to extract all the species names into a vector. We subset ...
# Blasthits$hits                 ... the list of hits, from which we choose ...
# Blasthits$hits[[i]]            ... the i-th hit, and get ...
# Blasthits$hits[[i]]$species    ... the species element from that.
# Subsetting FTW.

BLASTspecies <- character()
for (i in seq_along(BLASThits$hits)) {
    BLASTspecies[i] <-BLASThits$hits[[i]]$species
}

# You can confirm that BLASTspecies has the expected size.
length(BLASTspecies)

# Again, some species appear more than once, e.g. ...
sum(BLASTspecies == "Saccharomyces cerevisiae")

# ... corresponding to the five homologous gene sequences (paralogues) of yeast.

# Therefore we use unique() to throw out duplicates:
BLASTspecies <- unique(BLASTspecies)

length(BLASTspecies)
# i.e. we got rid of about two thirds of the hits.

# You should think about this: what is the biological interpretation of the
# finding that on average we have three sequences that are similar to Mbp1 in
# other species?


# =    4  Intersect GOLD and BLAST species  ====================================

# Now we can compare the two lists for species that appear in both sources: the
# simplest way is to use the set operation functions union(), intersection()
# etc. See here:
?union

MYSPEspecies <- intersect(GOLDspecies, BLASTspecies)

# Again: interpret this:
#  - what is the number of GOLDspecies?
#  - what is the number of BLAST species?
#  - how many species are present in both lists?
#  - what does it mean if a species is in GOLD but not in the BLAST list?
#  - what does it mean if a species has been found during BLAST, but it
#    is not in GOLD?


# =    5  Cleanup and finish  ==================================================

# One final thing: some of the species will be our so-called "reference" species
# which we use for model solutions and examples in the course. They are defined
# in the .utilities.R file of this project. We remove them from the list so that
# we don't inadvertently assign them.
#

REFspecies

MYSPEspecies <- sort(setdiff(MYSPEspecies, REFspecies))

# save(MYSPEspecies, file = "data/MYSPEspecies.RData")




# [END]
