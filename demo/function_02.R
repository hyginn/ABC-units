# tocID <- "demo/function_02.R"
#
#   Class session 2022-11-01
#   Starting code in class session 07, 2022-10-25
#
#  Boris Steipe (boris.steipe@utyoronto.ca)
#
#   Analysis motivated by Abou Chakra et al. (2021)
#    Control of tissue development and cell diversity by
#    cell cycle-dependent transcriptional filtering
#      https://elifesciences.org/articles/64951
#
# Version:   2.0
#
# Versions:
#            2.0 Complete rewrite 2022-11-02
#            1.0 Code protocol in class
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                          Line
#TOC> --------------------------------------------------------------
#TOC>   1        Context:                                         47
#TOC>   2        Supplementary data from paper                    65
#TOC>   3        Ensembl data:                                   173
#TOC>   3.1        Subsetting to remove irrelevant data          235
#TOC>   3.2        Subsetting to remove redundant data           253
#TOC>   3.3        Subsetting to remove unneeded columns         271
#TOC>   4        Milestone: storing the data                     298
#TOC>   5        Simple analysis                                 309
#TOC>   6        Distribution of annotations                     332
#TOC>   7        Distribution of gene lengths                    378
#TOC>   8        Finally: annotation of function                 514
#TOC>   8.1        From GOid to theme                            601
#TOC>   8.2        Wordclouds                                    645
#TOC> 
#TOC> ==========================================================================


# ==============================================================================



# =    1  Context:  ============================================================

# Last week we discussed the difficulties in annotating and differentiating the
# "function" of genes. Two genes can hardly be claimed to have the "same"
# function, but if all functions are unique, how could we say that a function is
# "enriched" in any set of genes"? Abou-Chakra et al. have taken an approach in
# which they group functions into "themes" and then see which of these themes
# could be observed more frequently in short and long genes. You should realize
# that the results of such an analysis depends crucially on how meaningful those
# themes are, and how well they group genes into meaningful categories that
# reflect similarities in their respective roles in the cell.

# This process is actually quite simple, and it gives us an opportunity to demo
# some fundamental data-science techniques. Exploring and validating the
# contents of data files, filtering and subsetting data, and matching elements
# from one dataset to another.


# =    2  Supplementary data from paper  =======================================

# Last week we ended off with looking at the semantics of two files we had
# downloaded as supplementary data from the Abou-Chakra et al. journal site.
# (For code to import these datasets, cf. function_01.R in this folder)

# dat1: contained keywords and themes...


# == Key data handling technique #1:  =================
#     - use   head()   to quickly inspect the contents of a large datafile
head(dat1)

#  Theme                  Keyword1
#  1        aging            adult
#  2     auditory         auditory
#  3     behavior         behavior   [...] (128 more columns)
#  4   blastocyst       blastocyst
#  5 blood-vessel     angiogenesis
#  [...]
#  (80 more rows)


# == Key data handling technique #2:  =================
#     - use   names()   early, to change column names to something sane

# dat2 mapped themes to actual GO terms in the "biological process" GO
# domain.
head(dat2)

# go_id                          pathway_name   pathway_name.1
# 1 GO:0055114    oxidation-reduction process        oxidation
# 2 GO:0009060            aerobic respiration             lung
# [...]
# (12,483 more rows)


# == Key data handling technique #3:  =================
#     - use   names()   EARLY, to change column names to something sane:

names(dat2) <- c("GOid", "GOdefinition", "Theme")

# It's especially important that columns with the same semantics don't have
# different names. That makes your code that much more readable.


# == Key data handling technique #4:  =================
#     - use   length(unique(...))   to count the number of unique items
#       in a dataset

# You can show that all dat2$GOid terms are unique ...
length(unique(dat2$GOid))                       # 12,485
length(unique(dat2$GOid)) == length(dat2$GOid)  # TRUE

# ... and that dat2 maps the 85 themes defined in dat1 to actual GOids.

length(unique(dat2$Theme))

# Let's have a brief look at those themes, to see if there are any problems with
# mapping one to the other. For brevity, we'll assign the relevant data items to
# two vectors: d1 and d2

d1 <- dat1$Theme            # all themes in dat1
d2 <- unique(dat2$Theme)    # all themes in dat2, each only once.


# == Key data handling technique #5:  =================
#     - the   %in%   operator ...

# Use the %in% operator to find all elements from one vector that are contained
# in the other. The output is a vector of logicals. Here we invert the result.
# Find all that are present in dat2$pathway_name.1 but NOT in dat1$Theme. Then
# we may subset dat2 with this vector to obtain the actual items:

sel <- ! (d2 %in% d1)
d2[sel]         # "growth-regulation"
                # "cell-adhesion"
                # "biological_process"
                # "muscle"

# Conversely: find all themes that were defined in dat 1, but not annotated to
# any of the GO_ids of dat2:

sel <- ! (d1 %in% d2)
d1[sel]         # "growth-regualtion"    <--! THIS IS A TYPO IN THE DATA!
                # "mammary-gland"
                # "muscle-control"
                # "pathogen"

# Fix the typo
which(dat1$Theme == "growth-regualtion")
dat1$Theme[which(dat1$Theme == "growth-regualtion")] <- "growth-regulation"

# After correcting the typo, that term is no longer different between the two
# files. While it is certainly unexpected that the "themes" present in those two
# files are not identical, we'll just remember that fact, in case we need it
# later.



# Now we have a better idea what to do next:
#   -  Navigate to biomart
#   -  Download a dataset of human genes, their gene lengths, their GO IDs ...
#   -  Add themes to their corresponding GO IDs
#   -  Check whether themes are different for short and for long genes.



# =    3  Ensembl data:  =======================================================

# - Accessed Ensembl as described  in the geneLengthFromBiomart.R script
#    http://uswest.ensembl.org/biomart/martview
# - Select human genes, annotate with:
#      GeneID and name,
#      Gene Start and end,
#      Gene type,
#      GO IDs, GO domain
#
# - Download 69,000 rows (130.9 MB)

# This file is surprisingly large for the 20,000 human genes we expect.

# I am not uploading this file to GitHub, but I will remove irrelevant data,
# redundant data and unneeded columns first. If you want to recreate the next
# three steps, you'll have to download the file from Ensebl. After the next
# three steps, you can read in the resulting (much smaller!) file from the
# course project.

# (In class, we originally worked with Transcript lengths - I now think that
# this is a mistake since transcripts are already spliced, and it is the
# pre-spliced version of the gene that is subject to transcriptional filtering.
# So I am re-doing this with gene-start and end. That's what the authors of
# the paper appear to have done, although the exact semantics of their data was
# not documented in the paper.)


inFile <- "~/Downloads/mart_export.txt"
hsGenes <- read.delim(inFile)

nrow(hsGenes)  # 1,594,546   Wow. Why?

names(hsGenes)
#  [1] "Gene.stable.ID"
#  [2] "Gene.start..bp."
#  [3] "Gene.end..bp."
#  [4] "Gene.name"
#  [5] "GO.domain"
#  [6] "GO.term.accession"
#  [6] "Gene.type"

# Clean up names
names(hsGenes) <- c("EnsembleID",
                    "start",
                    "end",
                    "name",
                    "GOdomain",
                    "GOid",
                    "type")

head(hsGenes)

#        EnsembleID start  end    name           GOdomain       GOid    type
# 1 ENSG00000210049   577  647   MT-TF molecular_function GO:0030533 Mt_tRNA
# 2 ENSG00000210049   577  647   MT-TF biological_process GO:0006412 Mt_tRNA
# 3 ENSG00000211459   648 1601 MT-RNR1 molecular_function GO:0003735 Mt_rRNA
# 4 ENSG00000211459   648 1601 MT-RNR1 cellular_component GO:0005840 Mt_rRNA
# 5 ENSG00000210077  1602 1670   MT-TV                               Mt_tRNA
# 6 ENSG00000210082  1671 3229 MT-RNR2 molecular_function GO:0003735 Mt_rRNA


# ==   3.1  Subsetting to remove irrelevant data  ==============================
#




# == Key data handling technique #6:  =================
#     - use ==, >, <, != etc. to subset a dataset to have only certain
#                             values in a column


# A: Only rows that are annotated to the "biological_process" GOdomain are
# relevant. All other rows can be removed.
hsGenes <- hsGenes[hsGenes$GOdomain == "biological_process" , ]
nrow(hsGenes)  # 538,219 ... was: 1,594,546



# ==   3.2  Subsetting to remove redundant data  ===============================



# == Key data handling technique #7:  =================
#     - use ! duplicated() to find unique items across more than one row ...

# B: remove duplicated combinations of gene name and GO id - these are
# redundant.
sel <- ! duplicated(hsGenes[ , c("name", "GOid") ])
hsGenes <- hsGenes[sel, ]

nrow(hsGenes)  # 164,436 ... was 538,219

# How many unique genes?
length(unique(hsGenes$name))  # 20,129  ... as we would expect


# ==   3.3  Subsetting to remove unneeded columns  =============================

# == Key data handling technique #8:  =================
#     - reassign a dataset with only a select set of columns to
#       reduce its size

# C: we don't actually ever use the Ensembl IDs, nor do we still need the
#    GOdomain column (all GOdomain values are now "biological process").
#    As well, the GOtermNames are already in dat2. So what we actually
#    need is only GeneName, start, end, GOid, and type. In that order.

hsGenes <- hsGenes[c("name", "start", "end", "GOid", "type")]

# Do all rows have GeneNames?

# == Key data handling technique #9:  =================
#     - use   sumt()   to count the number of TRUE elements in a logical vector.

sum(hsGenes$name == "")  # 513 rows do not have a gene name

# remove those:

sel <- hsGenes$name == ""
hsGenes <- hsGenes[! sel, ]

nrow(hsGenes)  # 163,923 ... was 164,436

# =    4  Milestone: storing the data  =========================================

# We store this data frame ...
# saveRDS(hsGenes, file = "demo/hsGenesGO.rds")

# ... and it's now only 1.2Mb in size. READ it in, and work along with the
# rest of the script!

hsGenes <- readRDS("demo/hsGenesGO.rds")


# =    5  Simple analysis  =====================================================


# == Key data handling technique #10:  =================
#     - perform some simple analysis to understand what you are dealing
#       with and can make decisions on how to proceed

# How many rows of data?
nrow(hsGenes) # 163,923

# How many unique genes?
length(unique(hsGenes$name)) # 20,128

# Remember: each row is a unique name / GOid combination ... but there
# may be several GOids annotated to each name.
# Average redundancy
nrow(hsGenes) / length(unique(hsGenes$name)) # 8.14

# On average we have a bit more than eight GOids annotated to each gene! But
# some may have only one, and some may have many more ... how are they
# distributed?


# =    6  Distribution of annotations  =========================================

# == Key data handling technique #11:  =================
#     - the many uses of   table()
#
#       table() is really one of the most valuable tools in the R toolchest.
#       Use it to:
#         - check whether there are typos in a column (those appear as unique
#           items);
#         - find the most frequent or least frequent items in data, ...
#         - ... or order data by frequency
#         - inspect the distribution of data item counts;
#         - feed a histogram or pie-chart;
#         - ...

# Distribution of number of GO ids annotated to each gene. Since every  name
# / GOid combination is unique, the number of GO ids annotated to each gene is
# obtained by inspecting how often a name appears in the data.
#

tGenes <- table(hsGenes$name)

# Let's see ...
summary(as.numeric(tGenes))  # Wow. Some gene has 213 annotations ...
#   Min.    1st Qu. Median  Mean    3rd Qu. Max.
#   1.000   2.000   5.000   8.144   9.000   213.000

# ... what is this gene?
tGenes[tGenes == 213]

# Ah. TNF. No wonder.

# Distribution?
boxplot(tGenes, col = "#c06c84")  # very skewed

hist(tGenes)

# The first 20 values count-by-count
hist(tGenes,
     breaks = 0:214,   # each bar is exactly one count different from the next
     xlim = c(0, 20),
     col = "#c06c84")




# =    7  Distribution of gene lengths  ========================================

# Distribution of gene lengths. We'll need that a lot, so we create a
# separate selection vector and keep it handy:

selUnique <- ! duplicated(hsGenes$name)


# == Key data handling technique #12:  =================
#     - validate.
#       Always validate!

sum(selUnique)                 # 20,128
length(unique(hsGenes$name))   # 20,128  MUST be the same


# === Histogram of lengths  (cf. Figure 8)



# == Key data handling technique #13:  =================
#     - add a column to a data frame by assigning to its not-yet-existing name.
#

hsGenes$length <- hsGenes$end - hsGenes$start

# Plot a histogram of the distribution:


# == Key data handling technique #14:  =================
#     - assign a   hist()    to a variable to work with the counts data later


myH <- hist(log10(hsGenes$length[selUnique]),
            col = "#fbcb7d",
            breaks = 40,
            xlim = c(0,8))

# Hm ... these are two distributions. Why? This NOT what Figure eight showed us.
# Have they made an additional, undocumented selection? What is in this peak? We
# can find out by inspecting the histogram, finding the boundaries of the bin we
# are curieous about, and then subsetting our data to that bin.

myH$counts[10:20]
#   [1]   23  268  516 1067  262  153   30  127   15  121  158

myH$breaks[10:20]
#   [1] 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7

# The peak of smaller genes is in bin 13 in the interval [2.0, 2.1]
sel <- log10(hsGenes$length) >= 2.0 &
       log10(hsGenes$length) <  2.1 &
       selUnique

# how many?
sum(sel)  # 1087

# what are these genes?
hsGenes$name[sel]  # sometimes you have to look at all of the data ...

# Fortunately Ensembl annotates gene types so we can get a high-level view:
unique(hsGenes$type[sel])

# [1] "snRNA"  "snoRNA" "miRNA"  "lncRNA" "scaRNA"

# RNA genes! It seems this peak contains RNA genes ... and it seems they have
# silently excluded all RNA genes without documenting it. Is this legitimate?
# Are RNA genes not subject to transcriptional filtering? After all, it is
# not "translational filtering" we are concerned with here? What do you think?

# Let's plot all RNA genes separately. First, we check what gene types we have
# in the first place:

(types <- unique(hsGenes$type))
#   [1] "Mt_tRNA"                          "protein_coding"
#   [3] "snRNA"                            "misc_RNA"
#   [5] "snoRNA"                           "lncRNA"
#   [7] "IG_V_gene"                        "IG_J_gene"
#   [9] "IG_C_gene"                        "IG_D_gene"
#   [11] "TR_V_gene"                        "TR_J_gene"
#   [13] "TR_D_gene"                        "scaRNA"
#   [15] "miRNA"                            "unprocessed_pseudogene"
#   [17] "ribozyme"                         "processed_pseudogene"
#   [19] "transcribed_processed_pseudogene" "transcribed_unitary_pseudogene"


# We can see three categories of types:

(typeRNA <-    c(types[grep("RNA", types)], "ribozyme"))
(typePseudo <- types[grep("pseudo", types)])
(typeProt   <- types[! (types %in% typeRNA | types %in% typePseudo)])

# (typeProt contains "normal" proteins, immunoglobulins and T-cell receptors.)

# Now replot our histogram:

selProtUnique <- ! duplicated(hsGenes$name) & (hsGenes$type %in% typeProt)
selRNAUnique  <- ! duplicated(hsGenes$name) & (hsGenes$type %in% typeRNA)

sum(selProtUnique)  # 17,550
sum(selRNAUnique)  #   2,571


# == Key data handling technique #15:  =================
#     - compare two distributions by overlaying their histograms
#

# First: define histogram-breaks that will cover BOTH distributions
myBreaks <- seq(0,8, length.out = 41)  # 40 intervals

# Second: plot your first distribution - the one in which you expect higher
# values to occur...
hist(log10(hsGenes$length[selProtUnique]),
     breaks = myBreaks,
     xlim = c(0, 8),
     col = "#fbcb7d")

# Third: plot your second distribution. Use the same breaks. Use the parameter
#        add = TRUE. And: use a transparent color.

hist(log10(hsGenes$length[selRNAUnique]),
     breaks = myBreaks,
     add = TRUE,
     col = "#d04c947f")  # the #......7f makes the color 50% transparent

# Well - there you have it. Figure eight is reproduced ... but not exactly. The
# authors have silently omitted RNA genes. Does this change their conclusions?
# Maybe. But unless you are capable of looking at the data yourself and making
# these observations, you can't even ask that question.

# The good thing is that R gives you all the tools you need to do this.

# But you actually have to do it. Yourself.



# =    8  Finally: annotation of function  =====================================

#
# Finally, we can turn back to our question of functional annotations.
#
# First: we subset our hsGenes to contain only gene types in typeProt ...

hsGenes <- hsGenes[hsGenes$type %in% typeProt, ]
nrow(hsGenes)   # 160,080 ... was: 163,923
length(unique(hsGenes$name))  # 17,550 ... was: 20,128

# Second: we make a version that contains each gene only once, because we need
# to figure out subsets of gene names according to gene length. That's easier if
# our selections don't always have to account for redundant genes.

sel <- ! duplicated(hsGenes$name)
hsUnique <- hsGenes[sel, ]
nrow(hsUnique)   # 17,550  as it should be

# Now we can order our table according to gene lengths!


# == Key data handling technique #16:  =================
#     -    order()    on a value to get a table sorted

# order() is one of those functions that's often neglected because you need a
# bit of mental effort to understand what it does. Don't be that way. order() is
# REALLY important.

ord <- order(hsUnique$length)

# ord is now a vector of the same length as hsUnique$length. And it contains the
# indices of the values of hsUnique$length IN ASCENDING ORDER. The short lengths
# come first:

head(ord)

# [1]  1303   580   475 14933 12276   579

# This means: the shortest gene is in row 1303, the second shortest in 580 etc.

tail(ord)

# [1]  8338  8806 16770  9830   917 14714

# This means: the longest gene is in row 14714, the second longest in 917 etc.

hsUnique[1303, ]    # shortest: probably an error in the database
hsUnique[14714, ]   # longest: 2.5 MB ... which is weird, because it (RBFOX1)
                    # is only 396aa long.

# But now we'll see why our order vector is so useful: it's trivially easy to
# use it to get the indices of the 5% shortest and longest genes.

n5 <- round(nrow(hsUnique) * 0.05)  # 878

l <- length(ord)
idxShort <- ord[         1:n5 ]  # indices of the 5% shortest genes ...
idxLong  <- ord[ (l-n5+ 1):l  ]  # ... and the 5% longest genes.

# Validate, validate, validate ...

# as boxplot ...
boxplot(log10(hsUnique$length[idxShort]),
        log10(hsUnique$length[ - c(idxShort, idxLong)]),
        log10(hsUnique$length[idxLong]),
        col = c("#BB44667F", "#fbcb7d", "#44AA887F"))
# Right? These should NOT overlap.

# as histogram ...
hist(log10(hsUnique$length[ - c(idxShort, idxLong)]),
     breaks = myBreaks,
     xlim = c(0, 8),
     col = "#fbcb7d")

hist(log10(hsUnique$length[idxShort]),
     breaks = myBreaks,
     add = TRUE,
     col = "#BB44667F")

hist(log10(hsUnique$length[idxLong]),
     breaks = myBreaks,
     add = TRUE,
     col = "#44AA887F")



# ==   8.1  From GOid to theme  ================================================

# For the functional annotations, we need to map our themes from dat2 to the GO
# ids we find in hsGenes. Actually: are all of the GOids even present in
# hsGenes?

sum( ! hsGenes$GOid %in% dat2$GOid) # 2,239 are not there! We'll have to be
                                    #   careful.



# == Key data handling technique #17:  =================
#     -    match()    cross-referecing values in two vectors
#                     match() is probably the second-most underappreciated
#                     function , after order(). Don't be that way. Learn
#                     to match! From the documentation:
#     match() returns a vector of the positions of (first) matches of its
#     first argument in its second argument.

# We want to find where the GOids of hsGene appear in dat2. Then we want to
# fetch the theme for that GOid, and put this into a separate column of
# hsGenes. Sound complicated? It isn't! match() does all the work.

hsGenes$themes <- dat2$Theme[ match(hsGenes$GOid, dat2$GOid) ]

# Validate ... just look at one:

i <- 77  # some random number
hsGenes[i, ]
dat2[dat2$GOid == hsGenes$GOid[i], ]

# Try a few more.


# Now we can pull out the themes for short genes ...

namesShort <- hsUnique$name[idxShort]
namesLong  <- hsUnique$name[idxLong]

themesShort <- hsGenes$themes[hsGenes$name %in% namesShort]
themesLong  <- hsGenes$themes[hsGenes$name %in% namesLong]



# ==   8.2  Wordclouds  ========================================================

# There are many ways to look at these collections of keywords. Many. Many are
# probably better than wordclouds. But we're doing wordclouds now. Because why
# not.

if (! requireNamespace("wordcloud", quietly = TRUE)) {
  install.packages("wordcloud")
}

tS <- table(themesShort)
tL <- table(themesLong)

wordcloud::wordcloud(names(tS),
                     tS^2,
                     max.words = 40,
                     random.order = FALSE,
                     colors = hcl.colors(40, "Temps"))
mtext("Wordcloud: GO themes in short Genes")


wordcloud::wordcloud(names(tL),
                     tL^2,
                     max.words = 30,
                     random.order = FALSE,
                     colors = hcl.colors(40, "Temps"))
mtext("Wordcloud: GO themes in long Genes")








# [END]
