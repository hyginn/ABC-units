# BIN-PPI-Analysis.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PPI-Analysis unit.
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




# ==============================================================================
#        PART FOUR: EXPLORE FUNCTIONAL EDGES IN THE HUMAN PROTEOME
# ==============================================================================

# In order for you to explore some real, biological networks, I give you a
# dataframe of functional relationships of human proteins that I have downloaded from the
# STRING database. The full table has 8.5 million records, here is a subset of
# records with combined confidence scores > 980

# The selected set of edges with a confidence of > 980 is a dataframe with about
# 50,000 edges and 6,500 unique proteins. You can load the saved dataframe here
# (and also read more about what the numbers mean at
# http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).

load("STRINGedges.RData")

head(STRINGedges)


# make a graph from this dataframe
?graph_from_data_frame

gSTR <- graph_from_data_frame(STRINGedges)

# CAUTION you DON'T want to plot a graph with 6,500 nodes and 50,000 edges -
# layout of such large graphs is possible, but requires specialized code. Google
# for <layout large graphs> if you are curious. Also, consider what one can
# really learn from plotting such a graph ...

# Of course simple computations on this graph are reasonably fast:

compSTR <- components(gSTR)
summary(compSTR) # our graph is fully connected!

dg <- degree(gSTR)
hist(log(dg), col="#FEE0AF")
# this actually does look rather scale-free

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#FEE0AF",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "6,500 nodes from the human functional interaction network")

# This looks very scale-free indeed.
#
# Now explore some more:

# === CLIQUES   ========
# Let's find the largest cliques. Remember: a clique is a fully connected
# subgraph, i.e. a subgraph in which every node is connected to every other.
# Biological complexes often appear as cliques in interaction graphs.

clique_num(gSTR)
# The largest clique has 63 members.

largest_cliques(gSTR)[[1]]

# Pick one of the proteins and find out what this fully connected cluster of 63
# proteins is (you can simply Google for the ID). Is this expected?



# === BETWEENNESS CENTRALITY   =======================================

# Let's find the nodes with the 10 - highest betweenness centralities.
#
BC <- centr_betw(gSTR)

# remember: BC$res contains the results
head(BC$res)

BC$res[1]   # betweeness centrality of node 1 in the graph ...
# ... which one is node 1?
V(gSTR)[1]

# to get the ten-highest nodes, we simply label the elements BC with their
# index ...
names(BC$res) <- as.character(1:length(BC$res))

# ... and then we sort:
sBC <- sort(BC$res, decreasing = TRUE)
head(sBC)

# This ordered vector means: node 3,862 has the highest betweeness centrality,
# node 1,720 has the second highest.. etc.

BCsel <- as.numeric(names(sBC)[1:10])
BCsel
# We can use the first ten labels to subset the nodes in gSTR and fetch the
# IDs...
ENSPsel <- names(V(gSTR)[BCsel])

# We are going to use these IDs to produce some output for you to print out and
# bring to class, so I need you to personalize ENSPsel with the following
# two lines of code:
set.seed(myStudentNumber)
ENSPsel <- sample(ENSPsel)

# We also need to remove the string "9606." from the ID:
ENSPsel <- gsub("9606\\.", "", ENSPsel)

# This is the final vector of IDs.:
ENSPsel

# Could you define in a short answer quiz what these IDs are? And what their
# biological significance is? I'm probably not going to ask you to that on the
# quiz, but I expect you to be able to.

#  Next, to find what these proteins are...

# We could now Google for all of these IDs to learn more about them. But really,
# googling for IDs one after the other, that would be lame. Let's instead use
# the very, very useful biomaRt package to translate these Ensemble IDs into
# gene symbols.

# == biomaRt =========================================================

# IDs are just labels, but for _bio_informatics we need to learn more about the
# biological function of the genes or proteins that we retrieve via graph data
# mining. biomaRt is the tool of choice. It's a package distributed by the
# bioconductor project. This here is not a biomaRt tutorial (that's for another
# day), simply a few lines of sample code to get you started on the specific use
# case of retrieving descriptions for ensembl protein IDs.

if (!require(biomaRt)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
  library("biomaRt")
}

# define which dataset to use ...
myMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# what filters are defined?
filters <- listFilters(myMart)
filters

# and what attributes can we filter for?
attributes <- listAttributes(myMart)
attributes

# Soooo many options - let's look for the correct name of filters that are
# useful for ENSP IDs ...
filters[grep("ENSP", filters$description), ]

# ... and the correct attribute names for gene symbols and descriptions ...
attributes[grep("symbol", attributes$description, ignore.case=TRUE), ]
attributes[grep("description", attributes$description, ignore.case=TRUE), ]


# ... so we can put this together: here is a syntax example:
getBM(filters = "ensembl_peptide_id",
      attributes = c("hgnc_symbol",
                     "wikigene_description",
                     "interpro_description",
                     "phenotype_description"),
      values = "ENSP00000000442",
      mart = myMart)

# A simple loop will now get us the information for our 10 most central genes
# from the human subset of STRING.

CPdefs <- list()  # Since we don't know how many matches one of our queries
# will return, we'll put the result dataframes into a list.

for (ID in ENSPsel) {
  CPdefs[[ID]] <- getBM(filters = "ensembl_peptide_id",
                        attributes = c("hgnc_symbol",
                                       "wikigene_description",
                                       "interpro_description",
                                       "phenotype_description"),
                        values = ID,
                        mart = myMart)
}

# So what are the proteins with the ten highest betweenness centralities?
#  ... are you surprised? (I am! Really.)

# Final task: Write a loop that will go through your list and
#    for each ID:
#    --  print the ID,
#    --  print the first row's symbol, and
#    --  print the first row's wikigene description.
#
# (Hint, you can structure your loop in the same way as the loop that
# created CPdefs. )

# Print the R code for your loop and its output for the ten genes onto a sheet
# of paper, write your student number and name on it, and bring this to class.


# = 1 Tasks




# [END]
