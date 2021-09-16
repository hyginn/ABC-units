# tocID <- "BIN-PPI-Analysis.R"
#
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-PPI-Analysis unit.
#
# Version:   1.4
#
# Date:     2017-08  -  2020-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.4    Update vector ID's for betweenness centrality.
#           1.3    Bugfix: called the wrong function on ENSPsel in l. 220
#           1.2    2020 Updates; Rewrite for new STRINg V11;
#                  Deprecate save()/load() for saveRDS()/readRDS()
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#                      use Biocmanager:: not biocLite()
#           1.0    First live version
#           0.1    First code copied from 2016 material.
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
#TOC>   Section  Title                                           Line
#TOC> ---------------------------------------------------------------
#TOC>   1        Setup and data                                    50
#TOC>   2        Functional Edges in the Human Proteome            86
#TOC>   2.1        Cliques                                        129
#TOC>   2.2        Communities                                    170
#TOC>   2.3        Betweenness Centrality                         184
#TOC>   3        biomaRt                                          231
#TOC>   4        Task for submission                              302
#TOC>
#TOC> ==========================================================================


# =    1  Setup and data  ======================================================


# Not surprisingly, the analysis of PPI networks needs iGraph:

if (! requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
# Package information:
#  library(help = igraph)       # basic information
#  browseVignettes("igraph")    # available vignettes
#  data(package = "igraph")     # available datasets

# In order for you to explore some real, biological networks, I give you a
# dataframe of functional relationships of human proteins that I have downloaded
# from the STRING database. The full table has 8.5 million records, here is a
# subset of records with combined confidence scores > 980

# The selected set of edges with a confidence of > 964 is a dataframe with about
# 50,000 edges and 8,400 unique proteins. Incidentaly, that's about the size of
# a fungal proteome. You can load the saved dataframe here (To read more about
# what the scores mean, see http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).

STRINGedges <- readRDS("./data/STRINGedges.rds")

head(STRINGedges)

# Note that STRING has appended the tax-ID for Homo sapiens - 9606 - to the
# Ensemble transcript identifiers that start with ENSP. We'll remove them:

STRINGedges$a <- gsub("^9606\\.", "", STRINGedges$a)
STRINGedges$b <- gsub("^9606\\.", "", STRINGedges$b)

head(STRINGedges)


# =    2  Functional Edges in the Human Proteome  ==============================


# There are many possibilities to explore interesting aspects of biological
# networks, we will keep with some very simple procedures here but you have
# to be aware that this is barely scratching the surface of possibilities.
# However, once the network exists in your computer, it is comparatively
# easy to find information online about the many, many options to analyze.


# Make a graph from this dataframe
?igraph::graph_from_data_frame

gSTR <- igraph::graph_from_data_frame(STRINGedges, directed = FALSE)

# CAUTION you DON'T want to plot a graph with 8,000 nodes and 50,000 edges -
# layout of such large graphs is possible, but requires specialized code. Google
# for <layout large graphs> if you are curious. Also, consider what one can
# really learn from plotting such a graph ...

# Of course simple computations on this graph are reasonably fast:

compSTR <- igraph::components(gSTR)
summary(compSTR) # our graph is fully connected!

hist(log(igraph::degree(gSTR)), col="#FEE0AF")
# this actually does look rather scale-free

(freqRank <- table(igraph::degree(gSTR)))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#FEE0AF",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "8,400 nodes from the human functional interaction network")

# This looks very scale-free indeed.

(regressionLine <- lm(log10(as.numeric(freqRank)) ~
                      log10(as.numeric(names(freqRank)) + 1)))
abline(regressionLine, col = "firebrick")

# Now explore some more:

# ==   2.1  Cliques  ===========================================================

# Let's find the largest cliques. Remember: a clique is a fully connected
# subgraph, i.e. a subgraph in which every node is connected to every other.
# Biological complexes often appear as cliques in interaction graphs.

igraph::clique_num(gSTR)
# The largest clique has 81 members.

(C <- igraph::largest_cliques(gSTR)[[1]])

# Pick one of the proteins and find out what this fully connected cluster of 81
# proteins is (you can simply Google for any of the IDs). Is this expected?

# Plot this ...
R <- igraph::induced_subgraph(gSTR, C) # a graph from a selected set of vertices

# color the vertices along a color spectrum
vCol <- rainbow(igraph::gorder(R)) # "order" of a graph == number of nodes

# color the edges to have the same color as the originating node
eCol <- character()
for (i in seq_along(vCol)) {
  eCol <- c(eCol, rep(vCol[i], igraph::gorder(R)))
}

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(R,
     layout = igraph::layout_in_circle(R),
     vertex.size = 3,
     vertex.color = vCol,
     edge.color = eCol,
     edge.width = 0.1,
     vertex.label = NA)
par(oPar)

# ... well: remember: a clique means every node is connected to every other
# node. We have 81 * 81 = 6,561 edges. This is what a matrix model of PPI
# networks looks like for large complexes.


# ==   2.2  Communities  =======================================================

set.seed(112358)                       # set RNG seed for repeatable randomness
gSTRclusters <- igraph::cluster_infomap(gSTR)
set.seed(NULL)                         # reset the RNG

igraph::modularity(gSTRclusters) # ... measures how separated the different
                                 # membership types are from each other
tMem <- table(igraph::membership(gSTRclusters))
length(tMem)  # About 700 communities identified
hist(tMem, breaks = 50, col = "skyblue")  # most clusters are small ...
range(tMem) # ... but one has > 200 members


# ==   2.3  Betweenness Centrality  ============================================

# Let's find the nodes with the 10 - highest betweenness centralities.
#
BC <- igraph::centr_betw(gSTR)

# remember: BC$res contains the results
head(BC$res)

BC$res[1]   # betweenness centrality of node 1 in the graph ...
# ... which one is node 1?
igraph::V(gSTR)[1]

# to get the ten-highest nodes, we simply label the elements of BC with their
# index ...
names(BC$res) <- as.character(1:length(BC$res))

# ... and then we sort:
sBC <- sort(BC$res, decreasing = TRUE)
head(sBC)

# This ordered vector means: node 3 has the highest betweenness centrality,
# node 721 has the second highest, etc.

(BCsel <- as.numeric(names(sBC)[1:10]))

# We can use the first ten labels to subset the nodes in gSTR and fetch the
# IDs...
(ENSPsel <- names(igraph::V(gSTR)[BCsel]))

# Task:
# =====
# IMPORTANT, IF YOU INTEND TO SUBMIT YOUR ANALYSIS FOR CREDIT
# We are going to use these IDs to produce some output for a submitted task:
# therefore I need you to execute the following line, note the "seal" that this
# returns, and not change myENSPsel later:

myENSPsel <- selectENSP(ENSPsel)

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

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
# Package information:
#  library(help = biomaRt)       # basic information
#  browseVignettes("biomaRt")    # available vignettes
#  data(package = "biomaRt")     # available datasets

# define which dataset to use ... this takes a while for download
myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

# what filters are defined?
( filters <- biomaRt::listFilters(myMart) )


# and what attributes can we filter for?
( attributes <- biomaRt::listAttributes(myMart) )


# Soooo many options - let's look for the correct name of filters that are
# useful for ENSP IDs ...
filters[grep("ENSP", filters$description), ]

# ... and the correct attribute names for gene symbols and descriptions ...
attributes[grep("symbol", attributes$description, ignore.case = TRUE), ]
attributes[grep("description", attributes$description, ignore.case = TRUE), ]


# ... so we can put this together: here is a syntax example:
biomaRt::getBM(filters = "ensembl_peptide_id",
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

for (ID in myENSPsel) {
  CPdefs[[ID]] <- biomaRt::getBM(filters = "ensembl_peptide_id",
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
# Write your thoughts about this group of genes.
#
# (Hint, you can structure your loop in the same way as the loop that
# created CPdefs. )

# Submit the "seal" for your ENSP vector, the ENSP vector itself, the R code
# for this loop and its output into your report if you are submitting
# anything for credit for this unit. Please read the requirements carefully.




# [END]
