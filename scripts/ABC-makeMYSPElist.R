# tocID <- "scripts/ABC-makeMYSPElist.R"
#
# Purpose:  Create a list of genome sequenced fungi with protein annotations and
#               Mbp1 homologues.
#
# Version: 1.5
#
# Date:    2016  09  -  2022  09
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions
#          1.5    2022 Maintenance
#          1.4    New retrieval logic
#          1.3    Rewrite to change datasource. NCBI has not been updated
#                   since 2012. Use ensembl fungi as initial source.
#          1.2    Change from require() to requireNamespace()
#          1.1.2  Moved BLAST.R to ./scripts directory
#          1.1    Update 2017
#          1.0    First code 2016
#
# TODO:
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
#TOC>   Section  Title                                    Line
#TOC> --------------------------------------------------------
#TOC>   1        The strategy                               55
#TOC>   2        PACKAGES AND INITIALIZATIONS               67
#TOC>   3        ENSEMBL FUNGI                              75
#TOC>   3.1        Import                                   78
#TOC>   4        BLAST SEARCH                              155
#TOC>   4.1        find homologous proteins                161
#TOC>   4.2        Identify species in "hits"              192
#TOC>   5        MERGE ENSEMBL AND BLAST RESULTS           282
#TOC>   6        STUDENT NUMBERS                           375
#TOC>
#TOC> ==========================================================================


# =    1  The strategy  ========================================================

# This script will create a list of "MYSPE" species and save it in an R object
# MYSPEspecies that is stored in the data subdirectory of this project from
# where it can be loaded. The strategy is as follows: we download a list of
# annotated fungal genomes from ensembl.fungi. All these are genome-sequenced
# species that have been annotated.
# Next we perform a BLAST search, to identify fungal species that have
# genes that are homologous to yeast MBP1.
#
# ...

# =    2  PACKAGES AND INITIALIZATIONS  ========================================

# httr provides interfaces to Webservers on the Internet
if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}


# =    3  ENSEMBL FUNGI  =======================================================


# ==   3.1  Import  ============================================================

# Navigate to https://fungi.ensembl.org and click on the link to the full
# list of all species: https://fungi.ensembl.org/species.html
# On the page, click on the spreadsheet symbol top right and choose
# "download whole table". The file will be named  "Species.csv", in your
# usual downloads folder. Move it to the data folder, and read it.

sDat <- read.csv("./data/Species.csv")
str(sDat)

# The most obvious way to partition these is according to Classification ...
# (poking around a bit in the UniProt taxonomy database shows that the
#  classification used here is the taxonomic rank of "order").
# how many classifications do we have?
length(unique(sDat$Classification))  # 66

# To have a good set for the class, we should have about 100.
# Let's see for which of these we can find Mbp1 homologues.
# First, we'll keep only the colums for name, classification, and taxID, and
# drop the rest ...
sDat <- sDat[ , c("Name", "Classification", "Taxon.ID")]
colnames(sDat) <- c("name", "order", "taxID")

# Next, we make an extra column: genus - the first part of the binomial name.
# We'll use the gsub() function, and for that we need a "regular expression"
# that matches to all characters from the first blank to the end of the string:
myPatt <- "\\s.*$"  # one whitespace (\\s) ...
                    # followed by any character (.) 0..n times (*) ...
                    # until the end of the string

# using gsub() we substitue all the matching characters with an empty
# string "" - this deletes the matching characters
# Test this:
gsub(myPatt, "", "Genus")                      # one word: unchanged
gsub(myPatt, "", "Genus species")              # two words: return only first
gsub(myPatt, "", "Genus species strain 123")   # many words: return only first

# apply this to the "name" column and add the result to the dataframe, as a
# separate column called "genus"
sDat$genus <- gsub(myPatt, "", sDat$name)

# what do we get?
c(head(unique(sDat$genus)), "...",
  tail(unique(sDat$genus)))  # inspect the first and last few genus.
                             # Note that there is a problem: we have an
                             # entry for "[Candida]" - the generic genus. There
                             # may be other special characters that will
                             # mess with our workflow. Let's drop those rows.
                             # (Always inspect your results!)
# what are they ?
sDat[ grepl("[^a-zA-Z]", sDat$genus) , ]

# drop them
sDat <- sDat[ ! grepl("[^a-zA-Z]", sDat$genus) , ]

length(table(sDat$genus))    # how many genus?
hist(table(sDat$genus), col = "#E9F4FF")      # Distribution ...
                                              # most genus have very few, but
                                              # some have very many species.
sort(table(sDat$genus), decreasing = TRUE)[1:10]  # Top ten...

# We should have at least one species from each taxonomic order, but we can
# add a few genus until we have about 100 validated species.

# Let's add a column for species, by changing our regular expression a bit,
# using ^ (start of string), \\S (NOT a whitespace),
# and + (one or more matches), capturing the match (...), and returning
# it as the substitution (\\1) ...

myPatt <- "^(\\S+\\s\\S+)\\s.*$"
sDat$species <- gsub(myPatt, "\\1", sDat$name)
head(sDat$species)

# And we reorder the columns, just for aesthetics:
sDat <- sDat[ , c("name", "species", "genus", "order", "taxID")]

# Final check:
any(grepl("[^a-zA-Z -.]", sDat$species)) # FALSE means no special characters

#
# Now we check which of these have Mbp1 homologues ...

# =    4  BLAST SEARCH  ========================================================


# We run a BLAST search to find all proteins related to yeast Mbp1 in any
# fungus. With the results, we'll annotate our sDat table.

# ==   4.1  find homologous proteins  ==========================================
#
# Use BLAST to fetch proteins related to Mbp1 and identify the species that
# contain them.

# Scripting against NCBI APIs is not exactly enjoyable - there is usually a fair
# amount of error handling involved that is not supported by the API in a
# principled way but requires rather ad hoc solutions. The code I threw together
# to make a BLAST interface (demo-quality, not research-quality) is in the file
# ./scripts/BLAST.R Feel encouraged to study how this works. It's a pretty
# standard task of communicating with servers and parsing responses - everyday
# fare in the bioinformatics lab. Surprisingly, there seems to be no good BLAST
# parser in currently available packages.
#
# DON'T use this for BLAST searches unless you have read the NCBI policy
# for automated tasks. If you indiscriminately pound on the NCBI's BLAST
# server, they will blacklist your IP-address. See:
# https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
#
# Use BLAST() to find yeast Mbp1 homologues in other fungi in refseq
# BLASThits <- BLAST("NP_010227",                  # Yeast Mbp1 RefSeq ID
#                    db = "refseq_protein",        # database to search in
#                    nHits = 3000,                 # 945 hits in 2020
#                    E = 0.01,                     #
#                    limits = "txid4751[ORGN]")    # = fungi
# saveRDS(BLASThits, file="data/BLASThits.2022-09-22.rds")
#
# NO NEED TO ACTUALLY RUN THIS: this was last run on 2022-09-22 and you can
#                               load the results from the data directory:
#
BLASThits <- readRDS(file = "data/BLASThits.2022-09-22.rds")

# ==   4.2  Identify species in "hits"  ========================================

# This is a very big list that can't be usefully analyzed manually. Here
# we are only interested in the species names that it contains.

# How many hits in the list?
length(BLASThits$hits)      # 1,262

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
    BLASTspecies[i] <- BLASThits$hits[[i]]$species
}

# You can confirm that BLASTspecies has the expected size.
length(BLASTspecies)

# if we delete some of these later on, we still want to remember which hit
# they came from. Thus we name() the elements with their index, which is the
# same as the index of the hit in BLASThits
names(BLASTspecies) <- 1:length(BLASTspecies)


# let's plot the distribution of E-values
eVals <- numeric()
for (i in seq_along(BLASThits$hits)) {
  eVals[i] <- BLASThits$hits[[i]]$E
}
range(eVals)
sum(eVals == 0)

# let's plot the log of all values > 0 to see how they are distributed
# plotting only one vecyor of numbers plots their index as x, and
# their value as y ...
plot(log(eVals[eVals > 0]), col = "#CC0000")

# This is very informative: I would suspect that the first fifty or so are
# virtually identical to the yeast protein, then we have about 900 hits with
# decreasing similarity, and then about 200 more that may actually be false
# positives. Also - we plotted them by index, that means the table is SORTED:
# Lower E-values strictly come before higher E-values.

# Again, some species appear more than once, e.g. ...
sum(BLASTspecies == "Saccharomyces cerevisiae")

# ... corresponding to the five homologous gene sequences (paralogues) of yeast.

# Therefore we remove duplicates. Removing duplicates will leave the FIRST in a
# list alone, and only remove the SUBSEQUENT ones. Remember, the list is sorted
# by E-values, thus removing duplicates retrieves the homologe from each
# species which is most similar to our query sequence, yeast Mbp1.
#
sel <- ! duplicated(BLASTspecies)  # positions of not-duplicated elements
head(sel)                           # all TRUE - first occurence
tail(sel)                           # all false - duplicates
BLASTspecies <- BLASTspecies[sel]   # filter the table and remove duplicates

length(BLASTspecies)
# i.e. we got rid of about two thirds of the hits.
tail(BLASTspecies)  # see how the names are useful!
                    # again - there are some special characters ...
                    # what are they?
BLASTspecies[grep("[^a-zA-Z ]", BLASTspecies)]

# remove the brackets ...
BLASTspecies <- gsub("\\[|\\]", "", BLASTspecies)
# drop any new duplicates ...
BLASTspecies <- BLASTspecies[ ! duplicated(BLASTspecies)]

# check the number again:
length(BLASTspecies)
# Think a bit about this: what may be the biological reason to find that
# on average, in 442 fungi across the entire phylogenetic tree, we have
# three sequences that are homologous to yeast Mbp1?

# Let's look at the distribution of E-values in this selection (Subsetting FTW):
# we plot all values that are TRUE in the vector "sel" that we created above,
# AND greater than 0
plot(log(eVals[sel & eVals > 0]), col = "#00CC00")



# =    5  MERGE ENSEMBL AND BLAST RESULTS  =====================================

# Next we add the BLAST result to our sDat dataframe. We'll store the index,
# the E-value, and the Query-bounds from which we can estimate which domains
# of Mbp1 are actually covered by the hit. (True orthologues MUST align with
# Mbp1's N-terminal APSES domain.)
#
# First we pull the hits we wanted from the BLASTspecies:
iHits <- as.numeric(names(BLASTspecies))
length(iHits)     # one index for each TRUE in sel

# add columns to sDat
l <- nrow(sDat)
sDat$iHit   <- numeric(l)  # index of the hit in the BLAST results
sDat$eVal   <- numeric(l)  # E-value of the hit
sDat$lAli   <- numeric(l)  # length of the aligned region

# extract and merge
for (iHit in iHits) {
  thisSp <- BLASThits$hits[[iHit]]$species
  sel <- sDat$species == thisSp

  sDat$iHit[sel]   <- iHit
  sDat$eVal[sel]   <- BLASThits$hits[[iHit]]$E
  sDat$lAli[sel]   <- BLASThits$hits[[iHit]]$lengthAli
}

# Are all reference species accounted for?
selA <- sDat$iHit != 0                 # all rows which matched to a BLAST hit
REFspecies %in% sDat$species[selA]     # yes, all there

selB <- sDat$species %in% REFspecies   # all rows which have one of REF species

sum(selA & selB)   # How many rows?

# sDat of course includes all duplicates. Some may be multiply sequenced, some
# may be different strains. We'll use the same strategy as before and keep
# only the best hit: order the rows by E-value, then drop all rows which
# are duplicated.


# drop all rows without BLAST hits ...
sDat <- sDat[ ! (sDat$iHit == 0) , ]

# order sDat by E-value ...
sDat <- sDat[order(sDat$eVal, decreasing = FALSE) , ]

# drop all rows with duplicated species ...
sDat <- sDat[ ! duplicated(sDat$species) , ]

# Lets look at the E-values ...
plot(log(sDat$eVal[sDat$eVal > 0]), col = "#0099AA")

# and alignment lengths ...
plot(sDat$lAli, col = "#00DDAA")

# How many ...
length(unique(sDat$name))
length(unique(sDat$species))
length(unique(sDat$genus))
length(unique(sDat$order))

# I need an extra species for admin purposes later on ...
sel <- grep("Sporothrix schenckii", sDat$species)
SPOSCdat <- sDat[sel, ]
sDat <- sDat[-sel, ]

# To get the final dataset, we remove the reference species with their
# entire orders ...
REForders <- unique(sDat$order[sDat$species %in% REFspecies])
sel <- sDat$order %in% REForders
REFdat <- sDat[sel , ]
sDat   <- sDat[ ! sel , ]

# REFdat should now contain only the ten REFspecies ...
( REFdat <- REFdat[REFdat$species %in% REFspecies , ] )

# ... but all of them
sum(REFspecies %in% REFdat$species)

# ... and we have enough left in sDat to prune sDat to unique genus
sDat <- sDat[ ! duplicated(sDat$genus) , ]
nrow(sDat)   # 99

# I add back "Sporothrix schenckii" ...
sDat <- rbind(SPOSCdat, sDat)

# ... and save for future use.
# saveRDS(sDat, file = "data/sDat.rds")
# saveRDS(REFdat, file = "data/REFdat.rds")



# =    6  STUDENT NUMBERS  =====================================================
#
# An asymmetric function to retrieve a MYSPE species from a student number.
# Asymmetric means: the student number should retrieve one and only one
# species, but it should not be possible to infer the student number from
# the species.
#
# sDat <- readRDS(file = "data/sDat.rds")

students <- read.csv("data/Roster-MGY441H1 F LEC5101 20229_Bioinformatics.csv")
print(summary(students$Student.Number), digits = 12)
sN <- students$Student.Number
sN <- sN[! is.na(sN)]
sN <- as.character(sN)
any(sN == "1003141593")
sN <- c("1003141593", sN)  # Artificial SN - will map to  "Sporothrix schenckii"

set.seed(112358)
theseSpecies <- sDat[sample(1:nrow(sDat)), ]
all(sort(theseSpecies$name) == sort(sDat$name))
nrow((theseSpecies))
(iX <- grep("Sporothrix schenckii", theseSpecies$name))
theseSpecies <- rbind(theseSpecies[iX, ], theseSpecies[-iX, ])
rndMin <-  992000000
rndMax <- 1020000000
N <- 10000
keys <- as.character(sample(rndMin:rndMax, N + 1000))
keys <- keys[! (keys %in% sN)]
keys <- keys[1:N]
keys[1:length(sN)] <- sN

nRep <- floor(N/nrow(theseSpecies))
MYSPEdat <- theseSpecies
for(i in 1:nRep) {
  MYSPEdat <- rbind(MYSPEdat, theseSpecies)
}
MYSPEdat <- MYSPEdat[1:N, ]
for (i in 1:N) {
  rownames(MYSPEdat)[i] <- digest::digest(keys[i], algo = "md5")
}
set.seed(NULL)
MYSPEdat <- MYSPEdat[sample(1:N), ]

# saveRDS(MYSPEdat, file = "data/MYSPEdat.rds")

# === validate ... (Recreate this if laqtecomer species are needed)
x <- character()
for (n in sN) {
  sp <- getMYSPE(n)
  if (length(sp) != 1) {
    stop(print(as.character(n)))
  } else {
    x <- c(x, sp)
  }
}
x
length(x)
length(unique(x))
getMYSPE("1003141593")


# === species for late-comers
y <- unique(MYSPEdat$species)
print(y[!(y %in% x)])


# === validate
l <- length(sN)
sp <- character(l)
for(i in 1:l) {
  sp[i] <- getMYSPE(sN[i])
}
any(duplicated(sp))
length(unique(sp))
which(! sDat$species %in% sp)  # these can be assigned to late-comers

# Done.

# [END]
