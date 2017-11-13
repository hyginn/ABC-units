# BIN-FUNC_Semantic_similarity.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-FUNC_Semantic_similarity unit.
#
# Version:  1.0
#
# Date:     2017  11  12
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    New code.
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
#TOC>   Section  Title                                          Line
#TOC> --------------------------------------------------------------
#TOC>   1        Preparations: Packages, AnnotationDB, Setup      39
#TOC>   2        Fetch GO Annotations                             89
#TOC>   3        Semantic Similarities                            98
#TOC>   4        GO Term Enrichment in Gene Sets                 116
#TOC>
#TOC> ==========================================================================


# =    1  Preparations: Packages, AnnotationDB, Setup  =========================


# GOSim is an R-package in the Bioconductor project.
if (! require(GOSim, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("GOSim")
  library(GOSim)
}
# Package information:
#  library(help = GOSim)       # basic information
#  browseVignettes("GOSim")    # available vignettes
#  data(package = "GOSim")     # available datasets


# GOSim loads human annotations by default. We load yeast annotations instead...
if (!require(org.Sc.sgd.db, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("org.Sc.sgd.db")
  library(org.Sc.sgd.db)
}

# org.Sc.sgd.db is a Bioconductor annotation database curated by SGD. Such
# databases exist for all model organisms. It's a kind of a fancy data frame
# from which we can get annotations by rows (genes) with the keys() funtion ...
keys(org.Sc.sgd.db)[1500:1510]

# ... and the types of available annotations with the columns() function
columns(org.Sc.sgd.db)

# Note that one of the columns is "GO" ... and we load that into the
# datastructures used by GOSim:

# Choose GOterms to use
setEvidenceLevel(evidences="all",
                 organism=org.Sc.sgdORGANISM,
                 gomap=org.Sc.sgdGO)

# Use Biological Process ontology
setOntology("BP", loadIC=FALSE)

# confirm that we loaded the correct ontology
head(get("gomap", envir=GOSimEnv))



# =    2  Fetch GO Annotations  ================================================


# All keys being used here are yeast systematic names.

# Get one set of annotations
getGOInfo(c("YDL056W"))  # Mbp1


# =    3  Semantic Similarities  ===============================================


# Get semantic similarities between genes
?getGeneSim

# There are _many_ different metrics of term similarity implemented
# in this package.

                                                     # Mbp1 and...
getGeneSim("YDL056W", "YLR182W", similarity = "OA")  # Swi6 - MCB complex
getGeneSim("YDL056W", "YER111C", similarity = "OA")  # Swi4 - collaborators
getGeneSim("YDL056W", "YBR160W", similarity = "OA")  # Cdc28 - mediator
getGeneSim("YDL056W", "YGR108W", similarity = "OA")  # Clb1 - antagonist
getGeneSim("YDL056W", "YLR079W", similarity = "OA")  # Sic1 - antagonist
getGeneSim("YDL056W", "YJL130C", similarity = "OA")  # Pgk1 - Gluconeogenesis


# =    4  GO Term Enrichment in Gene Sets  =====================================


# Calculating GO term enrichment in gene sets is done with the topGO package.
if (! require(topGO, quietly=TRUE)) {
  if (! exists("biocLite")) {
    source("https://bioconductor.org/biocLite.R")
  }
  biocLite("topGO")
  library(topGO)
}
# Package information:
#  library(help = topGO)       # basic information
#  browseVignettes("topGO")    # available vignettes
#  data(package = "topGO")     # available datasets


# Let's define a gene set: GOterm enrichment for G1/S switch activators:
mySet <- c("YFR028C", # Cdc14
           "YDL056W", # Mbp1
           "YLR182W", # Swi6
           "YER111C", # Swi4
           "YOR083W", # Whi5
           "YBR160W", # Cdc28
           "YMR199W", # Cln1
           "YPL256C", # Cln2
           "YAL040C") # Cln3

allGenes <- keys(org.Sc.sgd.db)
allGenes <- allGenes[grep("^Y", allGenes)]  # This is the context against which
                                            # we define enrichment

myEnr <- GOenrichment(mySet, allGenes)

sort(myEnr$p.values)  # Any significantly enriched terms? All of these are ...

#Yes: most significantly enriched is GO:0071931. What is this?
getGOTerm("GO:0071931")  # ... makes sense.

(fullSet <- myEnr$genes$`GO:0071931`)  # What genes are annotated to this term?

intersect(mySet, fullSet) # These are in both sets
setdiff(mySet, fullSet)   # These mySet members are not annotated to that term
setdiff(fullSet, mySet)   # These are annotated to that term but not in mySet.
                          # ... that's the most interesting set. From a set of
                          # genes we have identified a function that they
                          # share, and that shared function has allowed us
                          # to identify

# What are these genes?
# Select annotations from the annotation database:
select(org.Sc.sgd.db,
       keys = setdiff(fullSet, mySet),
       columns = c("COMMON", "DESCRIPTION"))

# Note that these annotations are partially redundant to several different
# aliases of the same three genes.



# [END]
