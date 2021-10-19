# 2021-10-12_In-Class_exploration.R
#
#         =====  T H E   E V E N   B E T T E R   A M I N O   A C I D =====
#
# Code and comments for BCH441 in-class exploration, Tuesday, 2021-10-12
# Explorers:  Jocelyn Nurtanto, Yuzi Li, and  Jerry Gu
# Scribe:     boris.steipe@utoronto.ca
#
# ==============================================================================
#
# In our last session we explored some properties of amino acids and noted that
# we can arrange them in a scatter-plot according to some properties. But can
# we also arrange them according to generic properties, i.e. taking all
# published property scales into account? We will try to use all tables from
# the seqinr package.

# First we load the package - this makes all datasets immediately available and
# we don't have to load them one by one.

library(seqinr)

# Determine what datasets are available
#
# Using "find in topic" ... "amino acid"
data(aacost)
data(aaindex)
data(pK)

# We note that datasets may be sorted in different ways: for example
# alphabetically by one letter code (A, C, D, E, ...) or three-letter code (Ala,
# Arg, Asn, Asp, ...) - this means we need to ensure and validate that amino
# acids are sorted in the same way.

# Build a datastructure ...
# rows: amino acids
# columns: properties

# Are all lists in aaindex organized in the same way?

refNames <- names(aaindex[[1]]$I) # Take the rownames of the first list item
                                  # index as a reference list

# Loop over each list in aaindex
for (i in 1:length(aaindex)) {
#   get the I-vector
  x <- aaindex[[i]]$I
#   get the names
  x <- names(x)
#   compare with the names of our reference list
#   the == and != operators are vectorized. Applying them to two vectors
#   gives TRUE or FALSE for each pair of elements. any() or all() can be
#   applied to logical vectors to anylise them and return a soingle result.
#   if (...) conditions evaluate only a single value and will throw a warning if
#   there is more than one.

  if (any(x != refNames)) {
    # There was at least one not-equal pair - so: complain
    print(sprintf("Problem in list %d: names don't match", i))
  }
}

# If we get here without identifying problems, it means all pairs of
# rownames match throughout the aainfex list.


# Next: what is the cvorrect syntax to add one vector (the "I" vector of
# one of the list elements) to our dataframe?
aaData <- as.data.frame(aaindex[[1]]$I) # Make a dataframe from the first index
aaData[,2] <- aaindex[[2]]$I            # ... add the secondf index

str(aaData)  # Confirm: we now have a two-column dataframe

# Next: add the rest ...
for (i in 3:length(aaindex)) {
  #   get the I-vector and write it into our dataframe
  aaData[,i] <- aaindex[[i]]$I
}

# Sanity check
plot(aaData[,37], aaData[,544])  # plot two arbitray inices against each other

# Looks good.

# We finished building our data structure ... but let's add the aacost table
# aacost is ordered differently:
rownames(aaData)
aacost[ , 1]

# using order(), applied to aacost - ordering the column with column-name
# "aaa"
sel <- order(aacost[ , "aaa"])  # alphebetic ordering of three-letter codes
aacost[sel, "aaa"] # applying the order vector sorts the column

# Is this the same order as refNames?
refNames == aacost[sel, "aaa"]  # Yes!

# add the data from column "tot" (i.e. total metabolic cost) after the
# last column of aaData
aaData[ , length(aaindex) + 1] <- aacost[sel, "tot"]

# Done.
str(aaData)  # A dataframe with 20 rows and 545 columns

# To answer the question "Which amino acids are similar to each other?" we
# need to reduce this 545-dimensional dataset to fewer dimensions, otherwise
# we will succumb to the "Curse of Dimensionality":
#
#    "in high dimensional data, however, all objects appear
#     to be sparse and dissimilar in many ways..."
#                   https://en.wikipedia.org/wiki/Curse_of_dimensionality
#
# A classic way to do this is Principal Component Analysis (PCA) ...
# (Principal components analysis)
#
# PCA expects objects in columns, properties in rows. Therefore we need to
# transpose our dataset:

aaPCA <- prcomp(t(aaData))

# This creates an error, because some of our indicews contain NA values!
# Which indices are this?

# We create a vector "sel" for which we check whether any element in each
# column is NA, and write FALSE if we encounter an NA, TRUE otherwise. We can
# then use this vector to subset ourt dataframe.

sel <- logical()

for (i in 1:ncol(aaData)) {         # for each index
  if (any(is.na(aaData[,i]))) {     #   if there is any NA value ...
    sel <- c(sel, FALSE)            #     add a FALSE element to the vector
  } else {                          #   else
    sel <- c(sel, TRUE)             #     add a TRUE element
  }
}

# Done. sel now subsets only the NA-free columns
545 - sum(sel)                      # 13 columns excluded

# Do the PCA ... use the prcomp() function
aaPCA <- prcomp(t(aaData[ ,sel]))   # PCA of the transposed, selected data set

str(aaPCA)   # structure of the result

plot(aaPCA)                         # plot the contributions of the
                                    # components to the variance

plot(aaPCA$rotation[ , 1],          # plot the first PC against the second PC
     aaPCA$rotation[ , 2],          # in a scatterplot, in an empty frame
     type ="n")                     # just to set up the coordinate system

text(aaPCA$rotation[ , 1],          # plot the names of the amino acids into
     aaPCA$rotation[ , 2],          # their respective (PC1, PC2) positions
     labels = rownames(aaPCA$rotation))

# PCA results are sensitive to the absolute numeric value of the features that
# we are comparing. The prcomp() function has an option scale. = TRUE that
# scales each row of features so that the variance of the value is 1.0  This
# ensures that each feature is given approximately equal weight

aaPCA <- prcomp(t(aaData[ ,sel]), scale. = TRUE)

plot(aaPCA)

plot(aaPCA$rotation[ , 1],
     aaPCA$rotation[ , 2],
     type ="n")
text(aaPCA$rotation[ , 1],
     aaPCA$rotation[ , 2],
     labels = rownames(aaPCA$rotation))


# Next we try to identify what the PCs correspond to. We see whether there are
# specific features that are highly correlated with the PCs

# ==== Rotation 1 ===================
#

(PC1 <- aaPCA$rotation[ , 1])  # Assign PC1

# The function cor() calculates Pearson coefficients of correlation
cor(PC1, aaData[ , 37]) # e.g. correlate PC1 against index 37


# Iterate over all columns and calculate correlations
cors <- numeric()

for (i in 1:ncol(aaData)) {
  cors[i] <- cor(PC1, aaData[ , i])
}

summary(cors)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's
# -0.54072 -0.13703  0.05654  0.03729  0.21349  0.59589       13
#
#  The max correlation is ~0.6. That is not very high. Which ijndex is it?

which(cors == max(cors, na.rm = TRUE))

aaindex[[504]]   # Linker propensity ???

cor(PC1, aaindex[[504]]$I) # Did we get the right index?

# Plot this ...
plot(aaPCA$rotation[ , 1],
     aaindex[[504]]$I,
     type ="n")
text(aaPCA$rotation[ , 1],
     aaindex[[504]]$I,
     labels = rownames(aaPCA$rotation))

# This is essentially a random correlation but for Cysteine ...


# ==== Rotation 2 ===================
#
# same process
PC2 <- aaPCA$rotation[ , 2]

cors2 <- numeric()

for (i in 1:ncol(aaData)) {
  cors2[i] <- cor(PC2, aaData[ , i])
}

summary(cors2)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's
# -0.95214 -0.56067 -0.12817 -0.05787  0.43046  0.94346       13

# Here we have quite strong correlations

which(cors2 == max(cors2, na.rm = TRUE))

aaindex[[148]]

# this index itself is correlated with many other indices

cor(PC2, aaindex[[148]]$I)   # confirmn that we have the right index

# Plot this too...
plot(aaPCA$rotation[ , 2],
     aaindex[[148]]$I,
     type ="n")
text(aaPCA$rotation[ , 2],
     aaindex[[148]]$I,
     labels = rownames(aaPCA$rotation))

# This correlates well with hydrophobicity measures. In this case the
# PC is to a certain degree interpretable - but this is not always the case
# with PCA (see the example of the first PC).






# [END]
