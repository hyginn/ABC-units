# tocID <- "scripts/ABC-makeScCCnet.R"
#
# Create a subnetwork of high-confidence yeast genes with a "mitotic cell cycle"
# GOSlim annotation.
#
# Boris Steipe for ABC learning units
#
# Notes:
#
#      The large source- datafiles are NOT posted to github. If you want to
#      experiment with your own code, download them and place them into your
#      local  ./data  directory.
#
#      STRING data source:
#        Download page:
# https://string-db.org/cgi/download.pl?species_text=Saccharomyces+cerevisiae
#        Data: (20.1 mb)
# https://stringdb-static.org/download/protein.links.full.v11.0/4932.protein.links.full.v11.0.txt.gz
#
#      GOSlim data source: (Note: this has moved from GO to SGD)
#        Info page: https://www.yeastgenome.org/downloads
#        Info page: http://sgd-archive.yeastgenome.org/curation/literature/
#        Data: (3 mb)
# http://sgd-archive.yeastgenome.org/curation/literature/go_slim_mapping.tab
#
#
# Version:  1.2
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2020 Update. GO Slim Yeast mow at SGD
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
#           1.0    First code copied from 2016 material.
#
# TODO:
#
# ==============================================================================
# SRCDIR <- "./instructor"


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                           Line
#TOC> ---------------------------------------------------------------
#TOC>   1        INITIALIZE                                        58
#TOC>   2        STRING FUNCTIONAL INTERACTION DATA                66
#TOC>   3        GOSlim FUNCTIONAL ANNOTATIONS                     96
#TOC>   3.1        Intersect interactions and annotations         122
#TOC>   4        DEFINE THE CELL-CYCLE NETWORK                    128
#TOC>
#TOC> ==========================================================================


# =    1  INITIALIZE  ==========================================================

SRCDIR <- "./data"
if (! requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}


# =    2  STRING FUNCTIONAL INTERACTION DATA  ==================================

# Read STRING Data (needs to be downloaded from database, see URL in Notes)
# The .gz compressed version is 20MB, the uncompressed versioj is 110MB -
# really not necessary to uncompress since readr:: can read from compressed
# files, and does so automatically, based on the file extension.
( fn <- file.path(SRCDIR, "4932.protein.links.full.v11.0.txt.gz") )
STR <- readr::read_delim(fn, delim = " ")

# Subset only IDs and combined_score column
STR <- STR[ , c("protein1", "protein2", "combined_score")]

# head(STR)
# sum(STR$combined_score > 909)  # 100270 edges
# subset for 100,000 highest confidence edges
STR <- STR[(STR$combined_score > 909), ]
head(STR)

# IDs are formatted like 4932.YAL005C ... drop the "4932." prefix
STR$protein1 <- gsub("^4932\\.", "", STR$protein1)
STR$protein2 <- gsub("^4932\\.", "", STR$protein2)
head(STR)

# get a vector of gene names in this list
myIntxGenes <- unique(c(STR$protein1, STR$protein2))  # yeast systematic gene
                                                      # names
length(myIntxGenes)
sample(myIntxGenes, 10)  # choose 10 at random (sanity check)


# =    3  GOSlim FUNCTIONAL ANNOTATIONS  =======================================
#
# Read GOSlim data  (needs to be downloaded from database, see URL in Notes)
( fn <- file.path(SRCDIR, "go_slim_mapping.tab") )

Gsl <- readr::read_tsv(fn,
                       col_names = c("ID",
                                     "name",
                                     "SGDId",
                                     "Ontology",
                                     "termName",
                                     "termID",
                                     "status"))

head(Gsl)

# What cell cycle names does it contain?
myGslTermNames <- unique(Gsl$termName)  # 169 unique terms
myGslTermNames[grep("cycle", myGslTermNames)]
# [1] "regulation of cell cycle"  "mitotic cell cycle"  "meiotic cell cycle"

# Choose "mitotic cell cycle" as the GOslim term to subset with

scCCgenes <- unique(Gsl$ID[Gsl$termName == "mitotic cell cycle"])
length(scCCgenes)  # 324 genes annotated to that term

# ==   3.1  Intersect interactions and annotations  ============================

sum(scCCgenes %in% myIntxGenes)  # 307 of these have high-confidence
#                                # functional interactions


# =    4  DEFINE THE CELL-CYCLE NETWORK  =======================================
#
# Define scCCnet ... the S. Cervisiae Cell Cycle network
# Subset all rows for which BOTH genes are in the GOslim cell cycle set
#
scCCnet <- STR[(STR$protein1 %in% scCCgenes) &
               (STR$protein2 %in% scCCgenes), ]

# How many genes are there?
length(unique(c(scCCnet$protein1, scCCnet$protein2)))  #283

# Each edge is listed twice - now remove duplicates.

# Step 1: make a vector: sort two names so the first one is alphabetically
#         smaller than the second one. This brings the two names into a defined
#         order. Then concatenate them with a "." - the resulting string
#         is always the same, for any order. E.g. c("A", "B") gives "A.B"
#         and c("B", "A") also gives "A.B". This identifies duplicates.

x <- apply(cbind(scCCnet$protein1, scCCnet$protein2),
           1,
           FUN = function(x) { return(paste(sort(x), collapse = ".")) })
head(x) # "YAL016W.YGR040W" "YAL016W.YOR014W" "YAL016W.YDL188C" ... etc.

sum(duplicated(x))  # 1453

# Step 2: drop all rows that contain duplicates in x
scCCnet <- scCCnet[! duplicated(x), ]

# Confirm we didn't loose genes
length(unique(c(scCCnet$protein1, scCCnet$protein2)))  # 283, no change
nrow(scCCnet)
# Network has 283 nodes, 1453 edges

saveRDS(scCCnet, file = "./data/scCCnet.rds")

# scCCnet <- readRDS("./data/scCCnet.rds")   # <<<- use this to restore the
                                             #      object when needed


# [END]
