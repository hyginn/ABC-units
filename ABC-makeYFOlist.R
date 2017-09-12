# ABC_makeYFOlist.R
#
# Purpose:  Create a list of genome sequenced fungi with protein annotations and
#               Mbp1 homologues.
#
# Version: 1.1
#
# Date:    2016  09 - 2017 08
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.1    Update 2017
# V 1.0    First code 2016
#
# TODO:
#   actually rerun for 2017
#   type out workflow
#
# ==============================================================================

# DO NOT  source()  THIS FILE!

# This file is code I provide for your deeper understanding of a process and
# to provide you with useful sample code. It is not actually necessary for
# you to run this code, but I encourage you to read it carefully and discuss
# if there are parts you don't understand.

# Run the commands that interact with the NCBI servers only if you want to
# experiment with the code and/or parameters. I have commented out those
# parts. If you simply want to reproduce the process you can simply
# load() the respective intermediate results.


# ==============================================================================
#        CREATING A YFO LIST
# ==============================================================================

# This script will create a list of "YFO" species and save it in an R object
# YFOspecies that is stored in the data subdirectory of this project from where
# it can be loaded.


# ==== GOLD species ============================================================
#
#  Fetch and parse genome data from the NCBI genome project database

# === Initialize
if (!require(httr)) { # httr provides interfaces to Webservers on the Internet
    install.packages("httr")
    library(httr)
}

# The URL where the genome data can be downloaded
URL <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt"

# Read the data directly from the NCBI ftp server and put it into a dataframe.
# This will take about a minute.

# GOLDdata <- read.csv(URL,
#                      header = TRUE,
#                      sep = "\t",
#                      stringsAsFactors = FALSE)
# save(GOLDdata, file="data/GOLDdata.RData")
load(file="data/GOLDdata.RData")


# What columns does the table have, how is it structured?
str(GOLDdata)

# What groups of organisms are in the table? How many of each?
table(GOLDdata$Group)

# What subgroups of fungi do we have?
table(GOLDdata$SubGroup[GOLDdata$Group == "Fungi"])

# How many of the fungi have protein annotations?
sum(GOLDdata$Proteins[GOLDdata$Group == "Fungi"] != "-")

# Get a subset of the data, with fungi that have protein annotations
GOLDfungi <- GOLDdata[GOLDdata$Group == "Fungi" &
                          GOLDdata$Proteins != "-" , ]

# check what we have in the table
head(GOLDfungi)

# For our purpose of defining species, we pick only the first two words from the
# "X.Organism.Name" column ... here is a function to do this:
#

makeBinomial <- function(s) {
    # input:
    #   s: a string which is expected to contain a binomial
    #      species name as the first two words, followed by other text
    # output:
    #   the first two words separeted by a single blank
    #
    x <- unlist(strsplit(s, "\\s+"))   # split second element on
    # one or more whitespace
    return(paste(x[1:2], collapse=" "))  # return first two elements
}

# iterate through GOLDdata and extract species names
GOLDspecies <- character()
for (i in 1:nrow(GOLDfungi)) {
    GOLDspecies[i] <- makeBinomial(GOLDfungi$X.Organism.Name[i])
}


# Species of great interest may may appear more than once, one for each sequenced strain: e.g. brewer's yeast:
sum(GOLDspecies == "Saccharomyces cerevisiae")

# Therefore we use the function unique() to throw out duplicates. Simple:
GOLDspecies <- unique(GOLDspecies)

length(GOLDspecies)
# i.e. we got rid of about half of the species.



# ==== BLAST species ===========================================================
#
# Use BLAST to fetch proteins related to Mbp1 and identifying the species that
# contain them.
#
# Scripting agains NCBI APIs is not exactly enjoyable - there is usually a fair
# amount of error handling involved that is not supported by the API in a
# principled way but requires rather ad hoc solutions. The code I threw
# together to make a BLAST interface for the course is in the file BLAST.R
# Feel encouraged to study how this works.

source("BLAST.R")   # load the function and its utilites

# Use BLAST() to find yeast Mbp1 homologues in other fungi in refseq
# hits <- BLAST("NP_010227",                  # Yeast Mbp1 RefSeq ID
#               nHits = 1000,                 # 633 hits in 2016
#               E = 0.01,                     #
#               limits = "txid4751[ORGN]")    # = fungi
# save(hits, file="data/BLASThits.RData")
load(file="data/BLASThits.RData")

# This is a very big list that can't be usefully analyzed manually. Here
# we are only interested in the species names that it contains.

# How many hits in the list?
length(hits$hits)

# Let's look at the first one
str(hits$hit[[1]])

# the species information is in the $def element - the definition line of the
# sequence record ... but we need a function to retrieve it. This one is a great
# example, because it is really messy. Have a look:
hits$hits[[1]]$def

# We get so many species because the exact same hit has been found by BLAST in a
# number of RefSeqs sequence records, all of which are different strains of
# yeast. We only need one of those though, any one, so we shall parse out the
# first one. For this, We will simply use the immensely versatile strsplit()
# function, split on square brackets, take the second element of the resulting
# array and run that through our makeBinomial function.

# define the function (i.e. execute the lines below)
parseDeflineSpecies <- function(s) {
    # input:
    #   s: a string which is expected to contain a binomial
    #      species name in square brackets embedded in other text
    # output:
    #   the species name found the first bracketed string
    #
    x <- unlist(strsplit(s, "\\]|\\["))  # split on "]" or "[" characters
    return(makeBinomial(x[2]))           # return Binomial name  from
                                         # second element of vector
}

#test it
parseDeflineSpecies(hits$hits[[1]]$def)
parseDeflineSpecies(hits$hits[[11]]$def)
parseDeflineSpecies(hits$hits[[111]]$def)

# now run a simple loop to extract all the species names into a vector
BLASTspecies <- character()
for (i in 1:length(hits$hits)) {
    BLASTspecies[i] <- parseDeflineSpecies(hits$hits[[i]]$def)
}

# You can confirm in the Values section of the Environment pane that
# BLASTspecies has the expected size. Again, species may appear more than once,
# e.g.
sum(BLASTspecies == "Saccharomyces cerevisiae")

# Therefore we use the function unique() to throw out duplicates. Simple:
BLASTspecies <- unique(BLASTspecies)

length(BLASTspecies)
# i.e. we got rid of about one third of the species.
#
# You should think about this: what does it mean that on average we have three
# hits by sequence similarity to Mbp1 in other species?


# ==== Intersecting BLAST and GOLD species lists ===============================


# Now we can compare the two lists for species that appear in both sources: the
# simplest way is to use the set operation functions union(), intersection()
# etc. See here:
?union

YFOspecies <- intersect(GOLDspecies, BLASTspecies)

# Just one final thing: some of the species will be our so-called "reference" species for which I will develop model solutions. I have defined them in the .utilities.R file to make them available for future purposes. separately and remove them from the list.
#

REFspecies

YFOspecies <- sort(setdiff(YFOspecies, REFspecies))

# save(YFOspecies, file = "data/YFOspecies.RData")


# [END]
