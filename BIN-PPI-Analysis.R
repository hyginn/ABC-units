# BIN-PPI-Analysis.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PPI-Analysis unit.
#
# Version:   1.0
#
# Date:     2017  08  - 2017 11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    First live version
#           0.1    First code copied from 2016 material.
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
#TOC>   Section  Title                                     Line
#TOC> ---------------------------------------------------------
#TOC>   1        Setup and data                              43
#TOC>   2        Functional Edges in the Human Proteome      80
#TOC>   2.1      Cliques                                    123
#TOC>   2.2      Communities                                164
#TOC>   2.3      Betweenness Centrality                     176
#TOC>   3        biomaRt                                    220
#TOC>   4        Task for submission                        291
#TOC> 
#TOC> ==========================================================================


# =    1  Setup and data  ======================================================


# Not surprisingly, the analysis of PPI networks needs iGraph:

if (!require(igraph, quietly=TRUE)) {
  install.packages("igraph")
  library(igraph)
}
# Package information:
#  library(help = igraph)       # basic information
#  browseVignettes("igraph")    # available vignettes
#  data(package = "igraph")     # available datasets

# In order for you to explore some real, biological networks, I give you a
# dataframe of functional relationships of human proteins that I have downloaded
# from the STRING database. The full table has 8.5 million records, here is a
# subset of records with combined confidence scores > 980

# The selected set of edges with a confidence of > 980 is a dataframe with about
# 50,000 edges and 6,500 unique proteins. Incidentaly, that's about the size of
# a fungal proteome. You can load the saved dataframe here (To read more about
# what the numbers mean, see http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).

load("./data/STRINGedges.RData")

head(STRINGedges)

# Note that STRING has appended the tax-ID for Homo sapiens - 9606 - to the
# Ensemble transcript identifiers that start with ENSP. We'll remove them:

STRINGedges$protein1 <- gsub("^9606\\.", "", STRINGedges$protein1)
STRINGedges$protein2 <- gsub("^9606\\.", "", STRINGedges$protein2)

head(STRINGedges)


# =    2  Functional Edges in the Human Proteome  ==============================


# There are many possibilities to explore interesting aspects of biological
# networks, we will keep with some very simple procedures here but you have
# to be aware that this is barely scratching the surface of possibilites.
# However, once the network exists in your computer, it is comparatively
# easy to find information nline about the many, many options to analyze.


# Make a graph from this dataframe
?graph_from_data_frame

gSTR <- graph_from_data_frame(STRINGedges, directed = FALSE)

# CAUTION you DON'T want to plot a graph with 6,500 nodes and 50,000 edges -
# layout of such large graphs is possible, but requires specialized code. Google
# for <layout large graphs> if you are curious. Also, consider what one can
# really learn from plotting such a graph ...

# Of course simple computations on this graph are reasonably fast:

compSTR <- components(gSTR)
summary(compSTR) # our graph is fully connected!

hist(log(degree(gSTR)), col="#FEE0AF")
# this actually does look rather scale-free

(freqRank <- table(degree(gSTR)))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#FEE0AF",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "6,500 nodes from the human functional interaction network")

# This looks very scale-free indeed.

(regressionLine <- lm(log10(as.numeric(freqRank)) ~
                      log10(as.numeric(names(freqRank)) + 1)))
abline(regressionLine, col = "firebrick")

# Now explore some more:

# ==   2.1  Cliques  ===========================================================

# Let's find the largest cliques. Remember: a clique is a fully connected
# subgraph, i.e. a subgraph in which every node is connected to every other.
# Biological complexes often appear as cliques in interaction graphs.

clique_num(gSTR)
# The largest clique has 63 members.

(C <- largest_cliques(gSTR)[[1]])

# Pick one of the proteins and find out what this fully connected cluster of 63
# proteins is (you can simply Google for any of the IDs). Is this expected?

# Plot this ...
R <- induced_subgraph(gSTR, C) # makes a graph from a selected set of vertices

# color the vertices along a color spectrum
vCol <- rainbow(gorder(R)) # gorder(): order of a graph = number of nodes

# color the edges to have the same color as the originating node
eCol <- character()
for (i in seq_along(vCol)) {
  eCol <- c(eCol, rep(vCol[i], gorder(R)))
}

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(R,
     layout = layout_in_circle(R),
     vertex.size = 3,
     vertex.color = vCol,
     edge.color = eCol,
     edge.width = 0.1,
     vertex.label = NA)
par(oPar)

# ... well: remember: a clique means every node is connected to every other
# node. We have 63 * 63 = 3,969 edges. This is what a matrix model of PPI
# networks looks like for large complexes.


# ==   2.2  Communities  =======================================================

set.seed(112358)
gSTRclusters <- cluster_infomap(gSTR)
modularity(gSTRclusters) # ... measures how separated the different membership
                         # types are from each other
tMem <- table(membership(gSTRclusters))
length(tMem)  # More than 2000 communities identified
hist(tMem, breaks = 50)  # most clusters are small ...
range(tMem) # ... but one has > 100 members


# ==   2.3  Betweenness Centrality  ============================================

# Let's find the nodes with the 10 - highest betweenness centralities.
#
BC <- centr_betw(gSTR)

# remember: BC$res contains the results
head(BC$res)

BC$res[1]   # betweeness centrality of node 1 in the graph ...
# ... which one is node 1?
V(gSTR)[1]

# to get the ten-highest nodes, we simply label the elements of BC with their
# index ...
names(BC$res) <- as.character(1:length(BC$res))

# ... and then we sort:
sBC <- sort(BC$res, decreasing = TRUE)
head(sBC)

# This ordered vector means: node 3,862 has the highest betweeness centrality,
# node 1,720 has the second highest.. etc.

(BCsel <- as.numeric(names(sBC)[1:10]))

# We can use the first ten labels to subset the nodes in gSTR and fetch the
# IDs...
(ENSPsel <- names(V(gSTR)[BCsel]))

# We are going to use these IDs to produce some output for a submitted task:
# so I need you to personalize ENSPsel with the following
# two lines of code:
set.seed(<myStudentNumber>) # enter your student number here
(ENSPsel <- sample(ENSPsel))

#  Next, to find what these proteins are...

# We could now Google for all of these IDs to learn more about them. But really,
# googling for IDs one after the other, that would be lame. Let's instead use
# the very, very useful biomaRt package to translate these Ensemble IDs into
# gene symbols.


# =    3  biomaRt  =============================================================


# IDs are just labels, but for _bio_informatics we need to learn more about the
# biological function of the genes or proteins that we retrieve via graph data
# mining. biomaRt is the tool of choice. It's a package distributed by the
# bioconductor project. This here is not a biomaRt tutorial (that's for another
# day), simply a few lines of sample code to get you started on the specific use
# case of retrieving descriptions for ensembl protein IDs.

if (!require(biomaRt, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("biomaRt")
  library(biomaRt)
}
# Package information:
#  library(help = biomaRt)       # basic information
#  browseVignettes("biomaRt")    # available vignettes
#  data(package = "biomaRt")     # available datasets

# define which dataset to use ...
myMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# what filters are defined?
(filters <- listFilters(myMart))


# and what attributes can we filter for?
(attributes <- listAttributes(myMart))


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


# =    4  Task for submission  =================================================


# Write a loop that will go through your personalized list of Ensemble IDs and
#    for each ID:
#    --  print the ID,
#    --  print the first row's HGNC symbol,
#    --  print the first row's wikigene description.
#    --  print the first row's phenotype.
#
# (Hint, you can structure your loop in the same way as the loop that
# created CPdefs. )

# Place the R code for this loop and its output into your report if you are
# submitting a report for this unit. Please read the requirements carefully.




# [END]
