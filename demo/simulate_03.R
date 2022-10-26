# simulate_03.R
#
#   Add randomization control (2022-10-26)
#   Next Phase - analysis - on 2022-10-18
#   (adding to Code, finished from class session 05, 2022-10-11)
#   (adding to Code from class session 04 - simulate_01)
#
#  Simulation motivated by Abou Chakra et al. (2021)
#    Control of tissue development and cell diversity by
#    cell cycle-dependent transcriptional filtering
#      https://elifesciences.org/articles/64951
#

# ==============================================================================

# Task:
#    I: Sketch the computational model of Abou-Chakra et al.
#   II: Write a simple implementation of an agent-based simulation
#  III: Experiment and analyse

#
# ==============================================================================
#    I: Sketch the computational model of Abou-Chakra et al
# ==============================================================================
#
#    Features of the model
#
# Globally: - define transcript rate (kb/hr)
#           - define target population
# Each cell: - is an agent
#            - has an ID
#            - has a lineage
#            - has genome
#            -   each gene
#              -   gene length
#              -   transcript count
#            - has cell cycle duration (hours)
#            - has a record of fractional cell cycle elapsed
#            - (no external factors)
#
# Each division:
#            - produces two new "active" cells
#              - defines cell cycle duration for each active cell
#              - store parent ID
#              - initializes transcript count (random apportioning)
#              - sets fractional cycle counter to zero
#            - sets the parent cell to "inactive"
#            - store child IDs
#
#
# Suitable data structures
#   - tree
#     - store children IDs
#
#
#  == Sketch simulation program flow ===========================================
#
#  Initialize first cell
#  For each active cell
#    update transcriptome
#    divide
#    break if number of cells is not less than target population
#  Update transcriptome one last time
#
#
#
# ==============================================================================
#   II: Write a simple implementation of an agent-based simulation
# ==============================================================================
#
# == Functions and datastructures. ====================
#
# == First phase: define what is required ======================================
#
#  Global parameters:
#   - the genome
#   - transcript rate
#
#  Datastructures
#   - myCells: the entire cell population
#   - cell: one row of a dataframe
#     - ID
#     - genes - transcript count
#     - cell cycle duration
#     - parent
#     - children
#       Note: each division creates two children from one parent. The
#             parent "ceases to exist". Only the children enter a next cycle.
#     - exists  (logical: whether this cell is currently alive)
#
#  Initialize the genome
#    - define each gene's ID and length
#
#  Initialize the first cell
#
#  Repeat
#    for each cell
#      if the size of the population has not been exceeded
#        update transcriptome
#        divide the cell
#        set parent to non-exist
#        increment cell counter by one
#      else: break
#
#
# == Second phase: datastructures and functions ================================
#
#  ==================
#  Global parameters:
#  ==================
MAXCELLS <- 10           # maximum number of cells existing in the experiment
TRATE <- 1               # - transcript rate (kb / h)
LCYC  <- 2               # - global (initial) cell-cycle length


# ===========
#  Functions:
# ===========

#  == newCells() ===============================================================
# Initialize the main data structure for the simulation.

newCells <- function(Nmax, lCyc, Genome = c(1, 2, 3)) {
  # Define the experiment's lineage

  cells <- list()

  # Define genome
  cells$genome <- Genome   # values are kb. Column numbers are Gene IDs.
  # The Cell table and the transcript table maintain a strict one-to-one
  # correspondence between rows. A cell ID is the row-number in which cell
  # information is stored in both the genome and the transcripts table.

  # Define transcript and cell tables. Note that the total number
  # of cells (existing + ancestral) is two times the number of existing cells.
  # But we need to add one extra in case the second child cell exceeds the
  # bounds.
  N <- (Nmax * 2) + 1

  # We keep the transcripts in a matrix. Each element [i, j] contains
  # the number of transcripts of gene j in cell i.
  cells$transcripts <- matrix(numeric(length(cells$genome) * N),
                              nrow = N)

  # The cell table holds both ancestral and existing cells.
  cells$cells <- data.frame(ID     = integer(N),  # unique ID
                            exists = logical(N),  # TRUE if cell exists
                            lCyc   = numeric(N),  # cycle length
                            parent = integer(N),  # ID of the parent cell
                            childA = integer(N),  # ID of first child cell
                            childB = integer(N))  # ID of second child cell

  #  Note: each division creates two children from one parent. The
  #        parent "ceases to exist". Only the children enter a next cycle.
  return(cells)
}

myCells <- newCells(MAXCELLS, LCYC)


#  == newID() ==================================================================
#  Return a new ID that is unique in the $cells$ID table.
newID <- function(cells) {
#' newID
#' @param   cells   The cells table. Must have a column called ID
#' @return int an integer ID that is one larger than the largest ID in the cells
#'                table.

  ID <- max(cells$cells$ID) + 1
  return(ID)
}

#  == initCell() ===============================================================
#  Add a first cell to $cells with initial values.
initCell <- function(cells) {
  # initialize a first cell
  cells$cells[1, ] <- data.frame(
    ID = newID(cells),          # unique ID
    exists = TRUE,              # TRUE if cell is  alive
    lCyc = LCYC,                # cycle length
    parent = 0,                 # ID of the parent cell
    childA = 0,                 # ID
    childB = 0)                 # ID
  return(cells)
}

#
#  Repeat
#    for each cell
#      if the size of the population has not been exceeded
#        update transcriptome
#
#  == updateTranscriptome() ====================================================
#  Bring a cell's transcriptome to the state at the end of the cell-cycle by
#  transcribing all of its genes for the duration of one cycle.
updateTranscriptome <- function(ID, cells) {
  lCyc <- cells$cells$lCyc[ID]                    # length of cycle
  lTrans <- lCyc * TRATE                          # length of transcript in
                                                  #   one cycle
  for (iGene in 1:length(cells$genome)) {         # for each gene
    nOld  <- cells$transcripts[ID, iGene]         # number of existing copies
    nNew  <- floor(lTrans / cells$genome[iGene])  # number of new copies
    cells$transcripts[ID, iGene] <- nOld + nNew   # update transcripts
  }
  return(cells)
}

#  == divideTranscriptome() ====================================================
#  Randomly apportion the existing transcriptome of a cell to its two child
#  cells.
divideTranscriptome <- function(IdP, IdA, IdB, cells) {
  for (iGene in 1:length(cells$genome)) {  # for each gene
    nP <- cells$transcripts[IdP, iGene]    # number of parent transcripts
    nA <- round(nP * runif(1))             # inheritance of child A
    nB <- nP - nA                          # inheritance of child A
    cells$transcripts[IdA, iGene] <- nA    # update transcripts
    cells$transcripts[IdB, iGene] <- nB    # update transcripts
  }
  return(cells)

}

#  == divideCell() =============================================================
# Divide the cell. Call divideTranscriptome().
divideCell <- function(IdP, cells) {
#          create two child cells and set their ID to IdA, IdB
#          also set parent IDs.

  # insert a first child cell into a new slot (row IdA)
  IdA <- newID(cells)   # Child A
  cells$cells[IdA, ] <- data.frame(
    ID=IdA, exists=TRUE, lCyc=LCYC, parent=IdP, childA=0, childB=0)

  # insert a second child cell into a new slot (row IdB)
  IdB <- newID(cells)   # Child B
  cells$cells[IdB, ] <- data.frame(
    ID=IdB, exists=TRUE, lCyc=LCYC, parent=IdP, childA=0, childB=0)

  # update parent's data
  cells$cells$exists[IdP] <- FALSE   # set parent to non-exist
  cells$cells$childA[IdP] <- IdA     # update child ID's
  cells$cells$childB[IdP] <- IdB

  # distribute parent's transcriptome to her children
  cells <- divideTranscriptome(IdP, IdA, IdB, cells)

  return(cells)
}


# == Third Phase: run the simulation and validate the code =====================

if (FALSE) {  # Do not execute this script part when the file is source()'d

  #  Initialize the cells
  myCells <- newCells(MAXCELLS, LCYC)

  #  Initialize the first cell
  myCells <- initCell(myCells)

  #  Repeat
  while (sum(myCells$cells$exists) <= MAXCELLS) {

    idx <- which(myCells$cells$exists)
    for (i in 1:length(idx)) {                     # for each existing cell
      # update transcriptome to prepare for division
      myCells <- updateTranscriptome(idx[i], myCells)
      # divide
      myCells <- divideCell(idx[i], myCells)
      if (sum(myCells$cells$exists) > MAXCELLS) {
        break
      }
    }
  }

  # Finish up: one last round of transcription

  idx <- which(myCells$cells$exists)
  for (i in 1:length(idx)) {                     # for each existing cell
    myCells <- updateTranscriptome(idx[i], myCells) # update transcriptome
  }

}

# Done.

# == Validate: =================================================================

# Increase length of cell cycle: longer transcripts appear. Correct
# Set cycle length to zero: no transcripts made. Correct.
# Set MAXCELLS to 0: only initial cell made. Correct.


LCYC <- 3
MAXCELLS <- 100

myCells <- newCells(MAXCELLS, LCYC)
myCells <- initCell(myCells)

while (sum(myCells$cells$exists) <= MAXCELLS) {

  idx <- which(myCells$cells$exists)
  for (i in 1:length(idx)) {                     # for each existing cell
    # update transcriptome to prepare for division
    myCells <- updateTranscriptome(idx[i], myCells)
    # divide
    myCells <- divideCell(idx[i], myCells)
    if (sum(myCells$cells$exists) > MAXCELLS) {
      break
    }
  }
}
idx <- which(myCells$cells$exists)
for (i in 1:length(idx)) {                     # for each existing cell
  myCells <- updateTranscriptome(idx[i], myCells) # update transcriptome
}

myCells$transcripts


# ==============================================================================
#    III: EXPERIMENT AND ANALYSE
# ==============================================================================

# ==  Simulate data for 5000 cells with a ten-gene genome:

LCYC <- 10
MAXCELLS <- 5000

myCells <- newCells(MAXCELLS, LCYC, Genome = 1:10)
myCells <- initCell(myCells)

while (sum(myCells$cells$exists) <= MAXCELLS) {
  idx <- which(myCells$cells$exists)
  for (i in 1:length(idx)) {
    myCells <- updateTranscriptome(idx[i], myCells)
    myCells <- divideCell(idx[i], myCells)
    if (sum(myCells$cells$exists) > MAXCELLS) { break }
  }
}
idx <- which(myCells$cells$exists)                # One last update
for (i in 1:length(idx)) {                        # for each existing cell
  myCells <- updateTranscriptome(idx[i], myCells) # update transcriptome
}

# Confirm
tail(myCells$transcripts)


# == Distributions

# Plot a boxplot()

# select all existing cells
sel <- which(myCells$cells$exists)

# When boxplot() is called on a matrix-like object, it creates one
# box-and-whisker object for each column.
boxplot(myCells$transcripts[sel, ])

# Compare:
summary(myCells$transcripts[sel, 1])

# Reproduce Abou Chakra et al. Figure 3A:

boxplot(list(c(myCells$transcripts[sel, 1:3]),   # short
             c(myCells$transcripts[sel, 4:7]),   # medium
             c(myCells$transcripts[sel, 8:10])), # long
        col = c("#EE984B", "#687AB2", "#B6949F"),
        bg = "#FAFAFA",
        outline = FALSE,
        xaxt = "n",
        whisklty = 1,
        boxcol = "#FFFFFF",
        border = "#AAAAAA",
        medlwd = 2,
        medcol=0)
legend("topright", legend = c("Short (1-3h)", "Medium (4-7h)", "Lomg (8-10h)"),
     fill = c("#EE984B", "#687AB2", "#B6949F"),
     bty = "n")

# Qualitatively, the author's observations are well reproduced.


# == Dimensionality Reduction ==================================================

# Abou Chakra et al. use the very capable Seurat package for an analysis of
# "cell diversity". Is it really true that transcriptional filtering alone can
# lead to groups of cells appearing, in a manner that suggests a process of
# differentiation? I would find that a very surprising result for such a
# stocahstic simulation, that has no strong inheritance within lineages.
#
# Our cells represent a 10-dimensional data-set, and in order to discover
# whether there is significant structure in this data, we can apply any of
# a number of dimensionality-reduction techniques, and / or clustering.
#
# Though Seurat:: makes a number of assumptions about the data that are not
# applicable to our simulation, the techniques themselves are easy enough
# to apply.
#
# Prepare data:

sel <- which(myCells$cells$exists)       # select "existing" cells
cellData <- myCells$transcripts[sel, ]   # subset data from transcripts table
cellData <- scale(cellData)              # scale columns to mean 0, sd 1
colMeans(cellData)                       # confirm: column means are now zero



# == Principal Components Analysis (PCA) =====================================

# PCA finds vectors on which one can "project" the data set, that give the
# projection a high variance.

pcDat <- prcomp(cellData)            # PCA
plot(pcDat)                          # plot variances of PCs ... this
                                     # is highly uncorrelated data. If
                                     # meaningful correlations could be found,
                                     # the left-most few Principal Components
                                     # would be much larger than the others.

plot(pcDat$x[ , c(1, 2)],            # Indeed, plotting the data on a scatter
     axes = FALSE, pch = 20,         # plot along the first two PCs shows that
     col=densCols(pcDat$x[,c(1,2)])) # no structure appears: there are no
                                     # separated clusters of points.

# You can experiment with plotting different scatterplots.
plot(pcDat$x[ , c(1, 3)],
     axes = FALSE, pch = 20,
     col = densCols(pcDat$x[ ,c(1,3)]))

plot(pcDat$x[ , c(2, 3)],
     axes = FALSE, pch = 20,
     col = densCols(pcDat$x[,c(2,3)]))


# == t-SNE (t-Stochastic Neighbour Embedding) ==================================
#
# t-SNE is an embedding method that preserves a kind of "distance" in high-
# dimensional space (a "t-statistic") as it embeds the data into a low-
# dimensional space. This means that data points that are close together in
# the original space remain close together after embedding. Just as with PCA
# the coordinates of the embedded result can not be directly interpreted -
# but the method is very useful, and widely used, to discover "structure" in
# the dataset - if there is any.
#
# cf. https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding

if (! requireNamespace("tsne", quietly = TRUE)) {
  install.packages("tsne")
}
library(tsne)

# t-SNE is quite slow for data sets as large as ours: all computational
# methods that take some form of "distance" between elements into account have
# (at least) an On^2 computational complexity. But we can run a valid
# experiment from a random subset of the data. Here we select just a thousand
# cells from our five thousand - which reduces the run-time by a factor of 25.

N <- 1000
rRow <- sample(1:nrow(cellData), N)   # random subset of row-indices of size N
smallDat <- cellData[rRow, ]          # subset the data

# The tsne function "iterates its embedding in cycles, and asks you to define
# a plotting function which it calls at each "epoch" of the embedding run.
# Here we just plot the data points, no axis labels, since the dimensions are
# arbitrary. The plotted coordinates are the two columns of the matrix which
# the algorithm returns.

# The tsne function "iterates its embedding in cycles, and you define a
# plotting function which it calls at each "epoch" of the embedding run. Here we
# just plot the data points, no axis labels, since the dimensions are arbitrary.
# The plotted coordinates are the two columns of the matrix which the algorithm
# returns.

myPlot <- function(x) {
  oPar <- par(mar=c(1, 1, 1, 1))     # reduce margin size
  plot(x[ , 1:2],                    # coordinates
       asp = 1,                      # use a square aspect ratio
       axes = FALSE,                 # no axes
       pch = 20,                     # round, filled dot
       col = densCols(x[ , 1:2]))    # cf. ?densCols
  par(oPar)                          # reset plot margins
}

# Try the function by plotting the first against the seventh column
# of the data set
myPlot(smallDat[ , c(1, 7)])

# Note again, that the data is essentially uncorrelated ...
cor(smallDat[ , 1], smallDat[ , 7])  # close to zero
# ... and discrete - and such discreteness obviously can be the reason
# why datasets may appear clustered.

# Running the t-SNE is simple;
tsneCells <- tsne::tsne(smallDat,
                       epoch_callback = myPlot, # Note: "myPlot", not
                                                #   "myPlot()" - we are
                                                #   passing a function, not
                                                #  its return value.
                       epoch = 10,              # number of iterations
                                                #   until a plot is made
                       perplexity = 20,         # size of the neighborhood
                                                #   the embedding considers
                       max_iter = 300)          # total number of
                                                #   iterations


# This is a stochastic process and every run will turn out somewhat differently,
# but you will probably find that while some apparent structure appears - as
# density maxima in the point clouds - in early iterations, after some one-
# hundred iterations or so, a much higher-quality solution is found that
# supports the idea that there actually is no inherent structure in this
# dataset. You can try this several times, and stop the process (by clicking on
# the little red STOP sign in the console menu bar) once there is no longer any
# significant change.

# Note that the points are not uniformly coloured! Given the discrete nature of
# the dataset, many of them overlap and appear darker.

# Are these "clusters"? Are they meaningful?

# What we would like to know is if the clusters that appear support the idea
# that there is a non-random process going on. Perhaps (and that would be
# exciting!) this simulation can show that transcriptional filtering can cause
# random fluctuations in early progenitors to get amplified to produce
# subpopulations of descendant cells? How could we pursue this question?

# One thing we could do is to trace back the lineages of cells in these clusters
# and see if they all derive from a common ancestor. Of course, they all derive
# from the single cell with which we started. But perhaps they share a common
# ancestor later on? Perhaps coloring the cells with 8 distinct colors
# corresponding to the fourth generation would show us such generational
# clustering. You can try that.

# But another way is to examine the properties of a random dataset.

# == A random control ...

# In order to examine this question, we could compare our simulated dataset to a random dataset. But if you think about this for a bit, you will realize that in our case, a random dataset is really hard to produce without fundamentally changing the nature of the experiment. We may get a result, but it may not be relevent. There is one thing we can do however - and this is a generic procedure that we can apply to very many different situations: shuffling our data.

    # =================================================================== #
    #                                                                     #
    #     Shuffling keeps all measures _on_ data points intact, and       #
    #     randomizes all associations _between_ data points.              #
    #                                                                     #
    # =================================================================== #

# Remember that forever. It is very useful.

# What association between data points would we want to remove? Well, we are
# thinking about lineages - and transcriptomes within lineages are partially
# inherited. A "transcriptome" here is the ten values within a single row. If we
# shuffle each column separately, there is no longer any association between a
# cell and its progenitor regarding how many transcripts were inherited. The
# transcriptomes are now completely random, althoug globally not a single
# transcript number has changed. They are just distributed differently. This is
# supe easy to code in R:


rndDat <- smallDat                          # copy the dataset

for (i in seq_len(ncol(rndDat))) {          # for each column
  rndDat[ , i] <- sample(rndDat[ , i])      # shuffle the values
}

# Done.
# Now we can run tsne with our shuffled data.

randCells <- tsne::tsne(rndDat,
                        epoch_callback = myPlot,
                        epoch = 10,
                        perplexity = 20,
                        max_iter = 300)

# Hm. What do you think?



# ==  Seurat complexity analysis ===============================================

# I think we are abusing Seurat when we apply it to this dataset - we might get
# it to produce plots, but the biological assumptions embedded into its
# workflows simply do not hold and that makes the results doubtful, IMO.

# That said, Seurat is an important analysis toolkit, especially for single-cell
# data and you are highly encouraged to work through the introductory tutorial
# at the Satija lab: https://satijalab.org/seurat/articles/get_started.html


# [END]
