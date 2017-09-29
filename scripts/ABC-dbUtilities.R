# ABC-dbUtilities.R

# database utilities for ABC learning units
#
# ==============================================================================
#


# ====== PACKAGES ==============================================================


if (! require("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
  library(jsonlite)
}


# ====== FUNCTIONS =============================================================


dbSanitizeSequence <- function(s, unambiguous = TRUE) {
  # Remove FASTA header lines, if any,
  # flatten any structure that s has,
  # remove all non-letters,
  # convert to uppercase.
  #
  # Parameters:
  #   s  chr  A DNA or protein sequence plus other characters
  #   unambiguous  bool  if TRUE, stop() if any letter remaining after
  #                      processing matches an ambiguity code. This is likely
  #                      due to inadvertently including meta-data, such as
  #                      a FASTA header, with the sequence.
  # Note: since U is an ambiguity code for amino acid sequences, you need
  #         to set unambiguous = FALSE to process RNA sequences with Uracil.
  # Value: chr   a valid, uppercase, amino acid sequence
  #

  s <- as.character(unlist(s))    # convert complex object to plain chr vector
  s <- unlist(strsplit(s, "\n"))  # split up at linebreaks, if any
  s <- s[! grepl("^>", s)]        # drop all lines beginning">" (FASTA header)
  s <- paste(s, collapse="")      # combine into single string
  s <- toupper(gsub("[^a-zA-Z]", "", s))
  if (unambiguous) {
    amb <- "([bjouxzBJOUXZ])"  # parentheses capture the match
    ambChar <- unlist(regmatches(s, regexec(amb, s)))[1]
    if (! is.na(ambChar)) {
      stop(paste("Input contains ambiguous codes(s): \"",
                 ambChar, "\".", sep=""))
    }
  }
  return(s)
}


dbConfirmUnique <- function(x) {
  # x is a vector of logicals.
  # returns x if x has exactly one TRUE element.
  # stop() otherwise.

  if (any(!is.logical(x))) {
    stop("PANIC: Input is not a boolean vector.")
  } else if (sum(x) == 0) {
    stop("PANIC: No match found.")
  } else  if (sum(x) > 1) {
    stop("PANIC: More than one match found.")
  } else {
    return(x)
  }
}


dbInit <- function() {
  # Return an empty instance of the protein database

  db <- list()

  db$protein <- data.frame(
    ID = numeric(),
    name = character(),
    RefSeqID = character(),
    UniProtID = character(),
    taxonomyID = numeric(),
    sequence = character(),
    stringsAsFactors = FALSE)

  db$taxonomy <- data.frame(
    ID = numeric(),
    species = character(),
    stringsAsFactors = FALSE)


  db$annotation <- data.frame(
    ID = numeric(),
    proteinID = numeric(),
    featureID = numeric(),
    start = numeric(),
    end = numeric(),
    stringsAsFactors = FALSE)

  db$feature <- data.frame(
    ID = numeric(),
    name = character(),
    description = character(),
    sourceDB = character(),
    accession = character(),
    stringsAsFactors = FALSE)

  return(db)
}


dbAutoincrement <- function(tb) {
  # Return a unique integer that can be used as a primary key
  # Value:
  #   num  a number one-larger than the largest current value in table$ID
  if (length(tb$ID) == 0) {
    return(1)
  } else {
    return(max(tb$ID) + 1)
  }
}


dbAddProtein <- function(db, jsonDF) {
  # Add one or more protein entries to the database db.
  # Parameters:
  #     db   list   a database created with dbInit()
  #     jsonDF  data frame  protein data imported into a data frame with
  #                           fromJSON()
  for (i in seq_len(nrow(jsonDF))) {
    x <- data.frame(ID          = dbAutoincrement(db$protein),
                    name        = jsonDF$name[i],
                    RefSeqID    = jsonDF$RefSeqID[i],
                    UniProtID   = jsonDF$UniProtID[i],
                    taxonomyID  = jsonDF$taxonomyID[i],
                    sequence    = dbSanitizeSequence(jsonDF$sequence[i]),
                    stringsAsFactors = FALSE)
    db$protein <- rbind(db$protein, x)
  }
  return(db)
}


dbAddFeature <- function(db, jsonDF) {
  # Add one or more feature entries to the database db.
  # Parameters:
  #     db   list   a database created with dbInit()
  #     jsonDF  data frame  feature data imported into a data frame with
  #                           fromJSON()
  for (i in seq_len(nrow(jsonDF))) {
    x <- data.frame(ID          = dbAutoincrement(db$feature),
                    name        = jsonDF$name[i],
                    description = jsonDF$description[i],
                    sourceDB    = jsonDF$sourceDB[i],
                    accession   = jsonDF$accession[i],
                    stringsAsFactors = FALSE)
    db$feature <- rbind(db$feature, x)
  }
  return(db)
}


dbAddTaxonomy <- function(db, jsonDF) {
  # Add one or more taxonomy entries to the database db.
  # Parameters:
  #     db      list         A database created with dbInit()
  #     jsonDF  data frame   Taxonomy data imported into a data frame with
  #                            fromJSON()
  for (i in seq_len(nrow(jsonDF))) {
    x <- data.frame(
      ID =  jsonDF$ID[i],
      species = jsonDF$species[i],
      stringsAsFactors = FALSE)
    db$taxonomy <- rbind(db$taxonomy, x)
  }
  return(db)
}

dbAddAnnotation <- function(db, jsonDF) {
  # Add one or more annotation entries to the database db.
  # Parameters:
  #     db   list   a database created with dbInit()
  #     jsonDF  data frame  annotation data imported into a data frame with
  #                           fromJSON()
  for (i in seq_len(nrow(jsonDF))) {

    sel <- jsonDF$pName[i] == db$protein$name
    sel <- dbConfirmUnique(sel)
    pID <- db$protein$ID[sel]

    sel <- jsonDF$fName[i] == db$feature$name
    sel <- dbConfirmUnique(sel)
    fID <- db$feature$ID[sel]

    x <- data.frame(ID        = dbAutoincrement(db$annotation),
                    proteinID = pID,
                    featureID = fID,
                    start     = as.integer(jsonDF$start[i]),
                    end       = as.integer(jsonDF$end[i]),
                    stringsAsFactors = FALSE)
    db$annotation <- rbind(db$annotation, x)
  }
  return(db)
}


# [END]
