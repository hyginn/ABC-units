# simulate_02.R
#
#
# Code from class session 05, 2022-10-11
#  (adding to Code from class session 04 - simulate_01)

# ==============================================================================

# Task:
#    I: Sketch the computational model of Abou-Chakra et al.
#   II: Write a simple implementation of an agent-based simulation
#

# ==============================================================================
#
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
# Suitable datastructures
#   - tree
#     - store children IDs
#
#  == Simulation program flow =====
#
#
#  Initialize first cell
#  For each active cell
#    divide
#
#    exit if number of cells is not less than target population
#
# == Functions and datastructures. ====================
#
# == First pass: define what is required
#
#  Global parameters:
#   - the genome
#
#  Datastructures
#   - cell: one row of a dataframe
#
#   - myCells: combines cell
#
#
#  Shared
#  Initialize genome
   myGenome <- c(0.5, 1.0)   # units: cell cycle time per transcript

#  Initialize first cell

  CELLID <- 0
  newCell <- function(thisID = (CELLID <- CELLID + 1),
                      thisG1 = 0,
                      thisG2 = 0) {
    cell <- data.frame(ID = thisID,
                       active = TRUE,
                       G1 = thisG1,
                       G2 = thisG2)
    return(cell)
  }

  myCells <- newCell()   # Initialize the population


#  For each active cell
#    divide
    # make it inactive
    # make two daughter cells
    CELLID <- CELLID + 2
#
#    exit if number of cells is not less than target population
#

myF <- function() {
  A <- 0
  x <- function() {
    y <- A
    A <<- A+1
    return(y)
  }
  return(x)
}

newID <- myF()

newID()

A <- -3
# ===========

newID <- (function() {
  A <- 0
  x <- function() {
    y <- A
    A <<- A+1
    return(y)
  }
  return(x)
})()



newID()






# [END]
