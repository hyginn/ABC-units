# tocID <- "./.utilities.R"
#
# Miscellaneous R code to suppport the project
#
# Version: 1.4
# Date:    2017-09 - 2020-09
# Author:  Boris Steipe
#
# V 1.4    Maintenance, and new validation utilities
# V 1.3.1  prefix Biostrings:: to subseq()
# V 1.3    load msa support functions
# V 1.2    update database utilities to support 2017 version of JSON sources
# V 1.1    2017 updates for ABC-units
# V 1.0    First code
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                       Line
#TOC> -----------------------------------------------------------
#TOC>   1        SCRIPTS TO SOURCE                             45
#TOC>   2        PACKAGES                                      51
#TOC>   3        SUPPORT FUNCTIONS                             62
#TOC>   3.1        objectInfo()                                65
#TOC>   3.2        biCode()                                    93
#TOC>   3.3        sameSpecies()                              127
#TOC>   3.4        pBar()                                     146
#TOC>   3.5        waitTimer()                                168
#TOC>   3.6        fetchMSAmotif()                            196
#TOC>   3.7        H() (Shannon entropy)                      240
#TOC>   4        DATA                                         254
#TOC>   4.1        REFspecies                                 256
#TOC>   5        FUNCTIONS TO CUSTOMIZE ASSIGNMENTS           271
#TOC>   5.1        getMYSPE()                                 274
#TOC>   5.2        selectPDBrep()                             283
#TOC>
#TOC> ==========================================================================


# =    1  SCRIPTS TO SOURCE  ===================================================

source("./scripts/ABC-dbUtilities.R")
source("./scripts/ABC-writeALN.R")
source("./scripts/ABC-writeMFA.R")

# =    2  PACKAGES  ============================================================

if (! requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest")
}

if (! requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}


# =    3  SUPPORT FUNCTIONS  ===================================================


# ==   3.1  objectInfo()  ======================================================
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


# ==   3.2  biCode()  ==========================================================
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


# ==   3.3  sameSpecies()  =====================================================
sameSpecies <- function(a, b) {
  # Parameters: a, b two vectors that contain
  # binomial species names and maybe additional strain information.
  # Value: a boolean vector, true where the species in a is the same as
  # the species in b.
  # Note: the usual vector recycling applies. Length is not checked.
  a <- gsub("^(\\S+\\s\\S+).*", "\\1", a)
  b <- gsub("^(\\S+\\s\\S+).*", "\\1", b)
  if (any(! grepl("^\\S+\\s\\S+$", a))) {
    stop("\"a\" contains elements that are not binomial names.")
  }
  if (any(! grepl("^\\S+\\s\\S+$", b))) {
    stop("\"b\" contains elements that are not binomial names.")
  }
  return(a == b)
}


# ==   3.4  pBar()  ============================================================
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


# ==   3.5  waitTimer()  =======================================================
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


# ==   3.6  fetchMSAmotif()  ===================================================
fetchMSAmotif <- function(ali, mot) {
  # Retrieve a subset from ali that spans the sequence in mot.
  # Biostrings package must be installed.
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

  return(Biostrings::subseq(ali@unmasked,
                start = which(x == motifStart)[1], # get the first position
                end   = which(x == motifEnd)[1]))
}


# ==   3.7  H() (Shannon entropy)  =============================================
H <- function(x, N) {
  # calculate the Shannon entropy of the vector x given N possible states
  # (in bits).
  # H(x) = - sum_i(P(x_i) * log2(P(x_i)); 0 * log(0) == 0
  t <- table(x)
  if (missing(N)) { N <- length(t) }
  if (length(t) > N ) { stop("N can't be smaller than observed states.") }
  h <- sum(- (t / length(x)) * log2(t / length(x)))
  return(h)
}



# =    4  DATA  ================================================================

# ==   4.1  REFspecies  ========================================================
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
                "Wallemia mellicola")


# =    5  FUNCTIONS TO CUSTOMIZE ASSIGNMENTS  ==================================

# ==   5.1  seal()  ========================================================
seal <- function(x.1L) { .Call(digest:::digest_impl,x.1L,3L,-1L,-0,-0,-0) }


# ==   5.1  getMYSPE()  ========================================================
getMYSPE <- function(x) {
  dat <- readRDS("./data/sDat.rds")
  map <- readRDS("./data/MYSPEmap.rds")
  key <- gsub(".+(....).$", "\\1", x)
  return(dat$species[map[key, "iMYSPE"]])
}


# ==   5.2  selectPDBrep()  ====================================================
selectPDBrep <- function(n, forCredit = FALSE) {
  # Select n PDB IDs from a list of high-resolution, non-homologous, single
  # domain, single chain structure files that represent a CATH topology
  # group.
  # Parameters:
  #   n     num     number of IDs to return
  #   seed  num     a seed for the RNG
  #
  # Value:          char  PDB IDs
  #
  # Note: the list is loaded from an .rds file in the "./data" directory.

  if (forCredit) {
    seed <- myStudentNumber
  } else {
    seed <- as.integer(Sys.time())
    cat("NOTE: This selection will not validate for a course submission.\n")
    cat("      If you intend to use it for an assignment task, invoke\n")
    cat("      this function like \"selectPDBrep(n, forCredit = TRUE)\".\n\n")
  }

  pdbRep <- readRDS("./data/pdbRep.rds")  # loads pdbRep

  if (n > length(pdbRep)) {
    stop(sprintf("There are only %d PDB IDs in the table to choose from.",
                 length(pdbRep)))
  }
  oldSeed <- .Random.seed
  set.seed(seed)
  PDBset <- sample(pdbRep, n)
  .Random.seed <- oldSeed
  return(PDBset)
}


# ==   5.2  selectChi2()  ====================================================
selectChi2 <- function() {
  # Select one random Amino acid from those that have a Chi2 angle

  oldSeed <- .Random.seed
  set.seed(myStudentNumber)
  AA <- sample(c("Asp", "Glu", "Phe", "His", "Ile", "Lys", "Leu",
                 "Met", "Asn", "Gln","Arg", "Trp", "Tyr"))
  .Random.seed <- oldSeed
  cat(sprintf("  Chi1/Ch2: Use \"%s\".  <%s>\n", AA[4], seal(AA)))
}


# ==   5.2  selectENSP()  ====================================================
selectENSP <- function(x) {
  oldSeed <- .Random.seed
  set.seed(myStudentNumber)
  x <- sample(x[order(x)])
  .Random.seed <- oldSeed
  cat(sprintf("seal: %s\n", seal(paste0(x,collapse=""))))
  return(x)
}



# [END]
