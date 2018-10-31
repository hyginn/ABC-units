# .utilities.R
#
# Miscellaneous R code to suppport the project
#
# Version: 1.3
# Date:    2017  09 - 2017  10
# Author:  Boris Steipe
#
# V 1.3    load msa support functions
# V 1.2    update database utilities to support 2017 version of JSON sources
# V 1.1    2017 updates for ABC-units
# V 1.0    First code
#
# ToDo:
# Notes:
#
# ==============================================================================

# ====== SCRIPTS =============================================================

source("./scripts/ABC-dbUtilities.R")
source("./scripts/ABC-writeALN.R")
source("./scripts/ABC-writeMFA.R")


# ====== SUPPORT FUNCTIONS =====================================================

objectInfo <- function(x) {
    # Function to combine various information items about R objects
    #
    # Input: an R object
    # Value: none - prints information as side-effect

    cat("object contents:")
    print(x, digits = 22)  # print value at maximal precision

    cat("\nstructure of object:\n")
    str(x)

    if (! is.list(x)) { # Don't use cat() if x is a list. cat() can't handle lists.
        cat("\nmode:   ", mode(x), "\n")
        cat("typeof: ", typeof(x), "\n")
        cat("class:  ", class(x), "\n")
    }

    # if the object has attributes, print them too
    if (! is.null(attributes(x))) {
        cat("\nattributes:\n")
        attributes(x)
    }
    # Done
}


biCode <- function(s) {
  # Make a 5 character "biCode" from a binomial name by concatening
  # the uppercased first three letter of the first word and the first
  # two letters of the second word. If there is only one word, we take the
  # first five characters from that. Outputs are padded with "." if necessary.
  # NAs in input are preserved.
  # Parameters:
  #    s   chr  vector of binomial species names
  # Value: chr  vector of biCodes, same length as s, NAs are preserved

  b <- character(length(s))
  s <- gsub("[^a-zA-Z ]", "", as.character(s)) # remove all non-alphabetic
                                               # characters except space
  s <- toupper(s)

  for (i in seq_along(s)) {
    x <- unlist(strsplit(s[i], "\\s+"))
    if (length(x) == 0) { # empty string
      x <- c("", "")
    } else if (length(x) == 1) { # only one string
      x <- c(substr(x, 1, 3), substr(x, 4, 5)) # 3 + 2 with whatever is there
    }
    x <- paste0(x[1:2], "...")  # pad strings

    b[i] <- paste0(substr(x[1], 1, 3), substr(x[2], 1, 2))
  }

  b[is.na(s)] <- NA  # recover NAs from input

  return(b)
}


pBar <- function(i, l, nCh = 50) {
  # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)[1]  # use only first, in case there are ties
      p1 <- paste(rep("#", p), collapse = "")
      p2 <- paste(rep("-", nCh - p), collapse = "")
      cat(sprintf("\r|%s%s|", p1, p2))
      flush.console()
    }
  }
  else { # done
    cat("\n")
  }
}


waitTimer <- function(t, nIntervals = 50) {
  # pause and wait for t seconds and display a progress bar as
  # you are waiting
  t <- as.numeric(t)

  if (t < 0.1) {return(invisible())}

  increment <- t / nIntervals

  bar <- "----:----|"  # One module for the progress bar:
  bar <- rep(bar, ceiling(nIntervals / 10))  # repeat,
  bar <- unlist(strsplit(bar, "")) # split into single characters,
  bar <- bar[1:nIntervals]  # truncate,
  bar <- paste(bar, collapse="") # and collapse.

  cat(sprintf("\nWaiting: |%s\n         |", bar))
  for (i in 1:(nIntervals - 1)) {
    Sys.sleep(increment)
    cat("=")
  }
  Sys.sleep(increment)
  cat("|\n\n")

  return(invisible())
}


fetchMSAmotif <- function(ali, mot) {
  # retrieve a subset from ali that spans the sequence in mot.
  # Parameters:
  #    ali        MsaAAMultipleAlignment object
  #    mot  chr   substring within ali
  # Value:  AAStringset

  if (class(ali) != "MsaAAMultipleAlignment" &&
      class(ali) != "MsaDNAMultipleAlignment" &&
      class(ali) != "MsaRNAMultipleAlignment") {
    stop("ali has to be an msa multiple alignment object.")
  }

  if (class(mot) != "character") {
    stop("mot has to be a character object.")
  }

  x <- gsub("-", "", as.character(ali))  # pure sequence, no hyphens

  idx <- grep(mot, x)[1] # first sequence containing mot. If no match,
                         # idx becomes NA
  if (is.na(idx)) {
    stop("mot is not a subsequence in ali.")
  }

  # Find the match range
  m <- regexpr(mot, x[idx])
  motifStart <- as.numeric(m)
  motifEnd <- attr(m, "match.length") + motifStart - 1

  # Count characters, skip hyphens ...
  x <- unlist(strsplit(as.character(ali)[idx], ""))
  x <- x != "-"
  x <- as.numeric(x)
  x <- cumsum(x)

  return(subseq(ali@unmasked,
                start = which(x == motifStart)[1], # get the first position
                end   = which(x == motifEnd)[1]))
}




# ====== PDB ID selection ======================================================

selectPDBrep <- function(n) {
  # Select n PDB IDs from a list of high-resolution, non-homologous, single
  # domain, single chain structure files that represent a CATH topology
  # group.
  # Parameters   n  num   number of IDs to return.
  # Value:          char  PDB IDs
  # Note: the list is loaded from an RData file in the data directory

  load("./data/pdbRep.RData")  # loads pdbRep
  if (n > length(pdbRep)) {
    stop(sprintf("You can select no more than %d IDs.", length(pdbRep)))
  }
  set.seed(as.numeric(Sys.time()))
  return(sample(pdbRep, n))
}


# ====== DATA ==================================================================


# 10 species of fungi for reference analysis.
# http://steipe.biochemistry.utoronto.ca/abc/index.php/Reference_species_for_fungi
REFspecies <- c("Aspergillus nidulans",
                "Bipolaris oryzae",
                "Coprinopsis cinerea",
                "Cryptococcus neoformans",
                "Neurospora crassa",
                "Puccinia graminis",
                "Saccharomyces cerevisiae",
                "Schizosaccharomyces pombe",
                "Ustilago maydis",
                "Wallemia mellicola"
)


# [END]
