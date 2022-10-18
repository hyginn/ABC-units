# simulate_02.R
#
#
# Code, finished from class session 05, 2022-10-11
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

newCells <- function(Nmax, lCyc) {
  # Define the experiment's lineage

  cells <- list()

  # Define genome
  cells$genome <- c(1, 2, 3)   # values are kb. Column numbers are Gene IDs.
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

# ...





# ==============================================================================
#    III: EXPERIMENT AND ANALYSE
# ==============================================================================


# ToDo: Analyze the change of transcriptome over the generations



# [END]
