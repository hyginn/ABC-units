# tocID <- "scripts/ABC-makeSTRINGedges.R"
#
# Create a subnetwork of high-confidence human STRING edges.
#
# Notes:
#
#      The large source- datafile is NOT posted to github. If you want to
#      experiment with the original data, download it and place it into your
#      local  ./data  directory.
#
#      STRING data source:
#        Download page:
# https://string-db.org/cgi/download.pl?species_text=Homo+sapiens
#        Data: (127.6 Mb)
# https://stringdb-static.org/download/protein.links.full.v11.0/9606.protein.links.full.v11.0.txt.gz
#
# Version:  1.0
#
# Date:     2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    Rewrite
#
# TODO:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                             Line
#TOC> -------------------------------------------------
#TOC>   1        Initialize                          44
#TOC>   2        Read STRING Data                    51
#TOC>   3        Define cutoff and subset            63
#TOC>   4        Drop  duplicates                   103
#TOC>   5        Simple statistics                  127
#TOC>   6        Write to file                      160
#TOC> 
#TOC> ==========================================================================


# =    1  Initialize  ==========================================================

if (! requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}


# =    2  Read STRING Data  ====================================================

# Read STRING Data (needs to be downloaded from database, see URL in Notes)
# The .gz compressed version is 127.6MB, the uncompressed version is probably
# 848 Mb. Fortunately readr:: can read from compressed
# files, and does so automatically, based on the file extension.
( fn <- file.path("~", "9606.protein.links.full.v11.0.txt.gz") )
STR <- readr::read_delim(fn, delim = " ")
nrow(STR)  #  11,759,454 rows
head(STR)


# =    3  Define cutoff and subset  ============================================

# approximate distribution of combined_score
hist(sample(STR$combined_score, 10000), breaks = 50, col = "#6699FF")

# Let's table the counts >= 850 and plot them for better resolution.

myTb <- table(STR$combined_score[STR$combined_score >= 850])
is.unsorted(as.integer(names(myTb)))  # Good - they are all in order

plot(myTb, type = "b", cex = 0.5, col = "#BB0000")
myTb[myTb == max(myTb)]  # Apparently there is an algorithmic effect that
                         # frequently assigns a combined score of 0.900

# Let's plot these counts as cumulative sums, in reverse order, scaled
# as combined scores.
myX <- 1 - (1:length(myTb)) / 1000   # x-values, decreasing
plot(myX,
     cumsum(myTb[length(myTb):1]),   # cumulative sum, decreasing
     xlim = c(1.0, 0.85),            # reverse x-axis
     type = "l",
     main = "STRING interactions for 9606 (top 600,000)",
     xlab = "combined_score",
     ylab = "cumulative counts",
     col = "#CC0000")
abline(h = seq(50000, sum(myTb), by = 50000), lwd = 0.5, col = "#DDDDFF")

# What's the cutoff for 100,000 edges?
which(cumsum(myTb[length(myTb):1]) >= 100000)[1] # p = 0.964

# confirm
sum(STR$combined_score >= 964) # 101,348
abline(v = 0.964, lwd = 0.5, col = "#DDDDFF")

# subset the table, and use only the protein IDs and the combined_score
STR <- STR[STR$combined_score >= 964,
            c("protein1", "protein2", "combined_score")]
colnames(STR) <- c("a", "b", "score")


# =    4  Drop  duplicates  ====================================================

# identify duplicate interactions by creating keys in a defined alphabetical
# sort order, then checking for  duplicated().
# e.g  if we have (X:U, U:X), we change U:X to X:U and now find that
# (X:U, X:U) has a duplicate.

AB <- STR$a < STR$b        # logical vector: genes we need to swap
tmp <- STR$b               # copy column b
STR$b[AB] <- STR$a[AB]     # copy a's into b
STR$a[AB] <- tmp[AB]       # copy tmp's into a
all(STR$a >= STR$b)        # confirm: TRUE

# now, make combined keys, like this:
paste0(STR$a[1:10], ":", STR$b[1:10])

tmp <- paste0(STR$a, ":", STR$b)
sum(duplicated(tmp)) # That's half of them ... i.e. STRING reports
                     # both a:b and b:a !

# drop all duplicated interactions from tmp
STR <- STR[ ! duplicated(tmp), ]   # 50,674 interactions remain


# =    5  Simple statistics  ===================================================

# how many unique genes?
length(unique(c(STR$a, STR$b)))   # 8,445

# how many self-edges?
sum(STR$a == STR$b)  # none

# log(rank) / log(frequency)
myTbl <- table(c(STR$a, STR$b))
myTbl <- myTbl[order(myTbl, decreasing = TRUE)]

hist(myTbl, breaks = 40, col = "#FFEEBB")

# number of singletons
sum(myTbl == 1) # almost a quarter

# maximum?
myTbl[which(myTbl == max(myTbl))]  # 9606.ENSP00000360532: 465
                                   # Google: CDC5L

# Zipf-plot
plot(log(1:length(myTbl)), log(as.numeric(myTbl)),
     type = "b", cex = 0.7,
     main = "STRINGedges - degrees",
     xlab = "log(rank)",
     ylab = "log(frequency)",
     col = "#FFBB88")

sprintf("Average number of interactions: %5.2f",
         nrow(STR) / length(unique(c(STR$a, STR$b))))


# =    6  Write to file  =======================================================

saveRDS(STR, file = "./data/STRINGedges.rds")

# STRINGedges <- readRDS("./data/STRINGedges.rds")  # use this to restore the
                                                    # object when needed


# [END]
