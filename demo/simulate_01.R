# simulate_01.R
#
# Code from class session 04, 2022-10-04

# ==============================================================================

# Task:
#    I: Sketch the computational model of Abou-Chakra et al.
#   II: Write a simple implementation of an agent-based simulation
#

# Simulate a coin flip ...

coin <- c("heads", "tails")
die <- 1:6
pFace <- c(0.1, 0.9) # really unfair coin

myAgent <- list(states = c("heads", "tails"),
                pState = c(0.1, 0.9))


x <- sample(myAgent$states, 3333, replace = TRUE, prob = myAgent$pState)
table(x)

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


  #  to be continued ...




# [END]
