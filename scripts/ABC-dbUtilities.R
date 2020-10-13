# tocID <- "scripts/ABC-dbUtilities.R"
#
# Purpose: Database utilities for ABC learning units.
#
# Version 2.2
#
# Date:     2017-11  -  2020-10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.2  Bugfixes
#           2.1  Add JSON export functions
#           2.0  Test all JSON import and prevent addition of duplicates. This
#                  is necessary for import of data from the public page
#           1.1  2020 Updates
#           1.0  Live version 2017
#
# Notes:
#   There are no functions to modify or delete entries. To do either,
#   recreate the database with correct data in the creation script. This is the
#   preferred way that ensures the entire database can be reproduced by
#   source()'ing its generating script.
#
#   Inserting data goes only through the very most minimal validation steps. For
#   production applications, more validation would need to be added, as well
#   as an overall validation of database integrity
#
# ToDo:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                   Line
#TOC> -------------------------------------------------------
#TOC>   1        INITIALISATIONS AND PARAMETERS            61
#TOC>   2        PACKAGES                                  66
#TOC>   3        FUNCTIONS                                 82
#TOC>   3.01       dbSanitizeSequence()                    85
#TOC>   3.02       dbConfirmUnique()                      120
#TOC>   3.03       dbInit()                               138
#TOC>   3.04       dbAutoincrement()                      178
#TOC>   3.05       dbAddProtein()                         191
#TOC>   3.06       dbAddFeature()                         227
#TOC>   3.07       dbAddTaxonomy()                        258
#TOC>   3.08       dbAddAnnotation()                      293
#TOC>   3.09       dbFetchUniProtSeq()                    340
#TOC>   3.10       dbFetchPrositeFeatures()               386
#TOC>   3.11       node2text()                            436
#TOC>   3.12       dbFetchNCBItaxData()                   448
#TOC>   3.13       UniProtIDmap()                         487
#TOC>   3.14       dbProt2JSON()                          526
#TOC>   3.15       dbSeq2JSON()                           611
#TOC>   3.16       dbRow2JSON()                           640
#TOC>   4        TESTS                                    660
#TOC> 
#TOC> ==========================================================================


# =    1  INITIALISATIONS AND PARAMETERS  ======================================

doTESTS <- FALSE  # run tests if TRUE


# =    2  PACKAGES  ============================================================


if (! requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}

if (! requireNamespace("xml2", quietly = TRUE)) {
  install.packages("xml2")
}


# =    3  FUNCTIONS  ===========================================================


# ==   3.01  dbSanitizeSequence()  =============================================
dbSanitizeSequence <- function(s, unambiguous = TRUE) {
  # Remove FASTA header lines, if any,
  # flatten any structure that s has,
  # remove all non-letters except "-" (gap) and "*" (stop),
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
  s <- toupper(gsub("[^a-zA-Z*-]", "", s))
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


# ==   3.02  dbConfirmUnique()  ================================================
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


# ==   3.03  dbInit()  =========================================================
dbInit <- function() {
  # Return an empty instance of the protein database
  # The schema is here:
  # https://docs.google.com/presentation/d/13vWaVcFpWEOGeSNhwmqugj2qTQuH1eZROgxWdHGEMr0

  db <- list()

  db$version <- "1.0"

  db$protein <- data.frame(
    ID = numeric(),
    name = character(),
    RefSeqID = character(),
    UniProtID = character(),
    taxonomyID = numeric(),
    sequence = character())

  db$taxonomy <- data.frame(
    ID = numeric(),
    species = character())

  db$annotation <- data.frame(
    ID = numeric(),
    proteinID = numeric(),
    featureID = numeric(),
    start = numeric(),
    end = numeric())

  db$feature <- data.frame(
    ID = numeric(),
    name = character(),
    description = character(),
    sourceDB = character(),
    accession = character())

  return(db)
}


# ==   3.04  dbAutoincrement()  ================================================
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


# ==   3.05  dbAddProtein()  ===================================================
dbAddProtein <- function(db, jsonDF) {
  # Add one or more protein entries to the database db if a protein with the
  # same name does not yet exist. This enforces that protein names are unique.
  # Parameters:
  #     db   list   a database created with dbInit()
  #     jsonDF  data frame  protein data imported into a data frame with
  #                           fromJSON()

  for (i in seq_along(jsonDF$name)) {
    isValid <- TRUE
    if (jsonDF$name[i] %in% db$protein$name) {
      cat(sprintf("Note: Protein No. %d in the input is \"%s\", but %s.\n",
                  i, jsonDF$name[i],
                  "a protein with this name already exists in the database. ",
                  "Skipping this input."))
      isValid <- FALSE
    }

    if (isValid) {
      if (length(jsonDF$name) == 1) { # jsonlite:: oversimplifies
        jsonDF$sequence <- paste(jsonDF$sequence, collapse = "")
      }
      x <- data.frame(ID          = dbAutoincrement(db$protein),
                      name        = jsonDF$name[i],
                      RefSeqID    = jsonDF$RefSeqID[i],
                      UniProtID   = jsonDF$UniProtID[i],
                      taxonomyID  = as.integer(jsonDF$taxonomyID[i]),
                      sequence    = dbSanitizeSequence(jsonDF$sequence[i]))
      db$protein <- rbind(db$protein, x)
    }
  }
  return(db)
}


# ==   3.06  dbAddFeature()  ===================================================
dbAddFeature <- function(db, jsonDF) {
  # Add one or more feature entries to the database db. Skip if a feature with
  # the same name already exists.
  # Parameters:
  #     db   list   a database created with dbInit()
  #     jsonDF  data frame  feature data imported into a data frame with
  #                           fromJSON()
  for (i in seq_along(jsonDF$name)) {
    isValid <- TRUE
    if (jsonDF$name[i] %in% db$feature$name) {
      cat(sprintf("Note: Feature No. %d in the input is \"%s\", but %s.\n",
                  i, jsonDF$name[i],
                  "a feature with this name already exists in the database. ",
                  "Skipping this input."))
      isValid <- FALSE
    }

    if (isValid) {
      x <- data.frame(ID          = dbAutoincrement(db$feature),
                      name        = jsonDF$name[i],
                      description = jsonDF$description[i],
                      sourceDB    = jsonDF$sourceDB[i],
                      accession   = jsonDF$accession[i])
      db$feature <- rbind(db$feature, x)
    }
  }
  return(db)
}


# ==   3.07  dbAddTaxonomy()  ==================================================
dbAddTaxonomy <- function(db, jsonDF) {
  # Add one or more taxonomy entries to the database db. Skip if species name
  # or taxonomy ID already exist in the database.
  # Parameters:
  #     db      list         A database created with dbInit()
  #     jsonDF  data frame   Taxonomy data imported into a data frame with
  #                            fromJSON()
  for (i in seq_along(jsonDF$species)) {
    isValid <- TRUE

    if (jsonDF$species[i] %in% db$taxonomy$species) {
      cat(sprintf("Note: Species No. %d in the input is \"%s\", but %s%s\n",
                  i, jsonDF$name[i],
                  "a species with this name already exists in the database. ",
                  "Skipping this input."))
      isValid <- FALSE
    } else if (jsonDF$ID[i] %in% db$taxonomy$ID) {
      cat(sprintf("Note: Taxonomy ID No. %d in the input is \"%d\", but %s%s\n",
                  i, jsonDF$ID[i],
                  "this taxonomy ID already exists in the database. ",
                  "Skipping this input."))
      isValid <- FALSE
    }
    if (isValid) {
      x <- data.frame(
        ID =  as.integer(jsonDF$ID[i]),
        species = jsonDF$species[i])
      db$taxonomy <- rbind(db$taxonomy, x)
    }
  }
  return(db)
}


# ==   3.08  dbAddAnnotation()  ================================================
dbAddAnnotation <- function(db, jsonDF) {
  # Add one or more annotation entries to the database db. Skip the entry if
  # it already exists in the database.
  # Parameters:
  #     db   list   a database created with dbInit()
  #     jsonDF  data frame  annotation data imported into a data frame with
  #                           fromJSON()
  for (i in seq_along(jsonDF$pName)) {
    isValid <- TRUE

    sel <- jsonDF$pName[i] == db$protein$name
    sel <- dbConfirmUnique(sel)    # Confirm that this protein ID exists
    pID <- db$protein$ID[sel]

    sel <- jsonDF$fName[i] == db$feature$name
    sel <- dbConfirmUnique(sel)    # Confirm that this feature ID exists
    fID <- db$feature$ID[sel]

    sel <- db$annotation$proteinID == pID &
           db$annotation$featureID == fID &
           db$annotation$start == as.integer(jsonDF$start[i]) &
           db$annotation$end   == as.integer(jsonDF$end[i])

    if (any(sel)) {
      cat(sprintf("Note: annotation No. %d in the input has %s%s\n",
                  i,
                  "the same protein name, feature name, start, and end ",
                  "as one that already exists in the database. ",
                  "Skipping this input."))

      isValid <- FALSE
    }

    if (isValid) {
      x <- data.frame(ID        = dbAutoincrement(db$annotation),
                      proteinID = pID,
                      featureID = fID,
                      start     = as.integer(jsonDF$start[i]),
                      end       = as.integer(jsonDF$end[i]))
      db$annotation <- rbind(db$annotation, x)
    }
  }
  return(db)
}


# ==   3.09  dbFetchUniProtSeq()  ==============================================
dbFetchUniProtSeq <- function(IDs) {
  # Fetch a protein sequence from UniProt.
  # Parameters:
  #     IDs   char  a vector of UniProt IDs (accession number)
  # Value:
  #     char        a vector of the same length as ID. It contains
  #                 sequences where the retrieval was successful, NA where
  #                 it was not successful. The elements are named with
  #                 the ID, the header lines are set as attribute "header"


  BASE <- "http://www.uniprot.org/uniprot/"

  sq <- character()
  hd <- character()
  for (i in seq_along(IDs)) {
    URL <- sprintf("%s%s.fasta", BASE, IDs[i])
    response <- httr::GET(URL)
    if (httr::status_code(response) == 200) {
      s <- as.character(response)
      s <- unlist(strsplit(s, "\n"))
      x <- dbSanitizeSequence(s)
    } else {
      s <- ""
      x <- NA
    }
    hd[i] <- s[1]
    sq[i] <- x
  }
  names(sq) <- IDs
  attr(sq, "headers") <- hd

  return(sq)
}

if (FALSE) {
  inp <- c("P79073", "P0000000", "A0A1W2TKZ7")
  s <- dbFetchUniProtSeq(inp)
  s[1:3]
  str(s)
  attr(s, "headers")[1]
}



# ==   3.10  dbFetchPrositeFeatures()  =========================================
dbFetchPrositeFeatures <- function(ID) {
  # Fetch feature annotations from ScanProsite.
  # Parameters:
  #     ID   char   a UniProt ID (accession number)
  # Value:
  #     data frame     uID    char  UniProt ID
  #                    start  num   start of motif
  #                    end    num   end of motif
  #                    psID   char  PROSITE motif ID
  #                    psName char  PROSITE motif name
  #                    psSeq  char  sequence annotated to the feature
  # If the operation is not successful, a 0-length data frame is returned.

  URL <- "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"

  response <- httr::POST(URL,
                         body = list(meta = "opt1",
                                     meta1_protein = "opt1",
                                     seq = ID,
                                     skip = "on",
                                     output = "tabular"))

  myFeatures <- data.frame()
  if (httr::status_code(response) == 200) {

    lines <- unlist(strsplit(httr::content(response, "text"), "\\n"))

    patt <- sprintf("\\|%s\\|", ID)
    lines <- lines[grep(patt, lines)]

    for (line in lines) {
      tokens <- unlist(strsplit(line, "\\t|\\|"))
      myFeatures <- rbind(myFeatures,
                          data.frame(uID   =  tokens[2],
                                     start =  as.numeric(tokens[4]),
                                     end   =  as.numeric(tokens[5]),
                                     psID  =  tokens[6],
                                     psName = tokens[7],
                                     psSeq  = tokens[11]))
    }
  }
  return(myFeatures)
}

if (FALSE) {
  dbFetchPrositeFeatures("P33520")  # RES1_SCHPO

}

# ==   3.11  node2text()  ======================================================
node2text <- function(doc, tag) {
  # an extractor function for the contents of elements
  # between given tags in an XML response.
  # Contents of all matching elements is returned in
  # a vector of strings.
  path <- paste0("//", tag)
  nodes <- xml2::xml_find_all(doc, path)
  return(xml2::xml_text(nodes))
}


# ==   3.12  dbFetchNCBItaxData()  =============================================
dbFetchNCBItaxData <- function(ID) {
  # Fetch feature taxID and Organism from the NCBI.
  # Parameters:
  #     ID   char   a RefSeq ID (accession number)
  # Value:
  #     data frame     taxID      num    NCBI taxID
  #                    organism   char   organism for this taxID
  # If the operation is not successful, a 0-length data frame is returned.

  eUtilsBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  URL <- paste(eUtilsBase,
               "esearch.fcgi?",
               "db=protein",
               "&term=", ID,
               sep="")
  myXML <- xml2::read_xml(URL)
  GID <- node2text(myXML, "Id")

  URL <- paste0(eUtilsBase,
                "esummary.fcgi?",
                "db=protein",
                "&id=",
                GID,
                "&version=2.0")
  myXML <- xml2::read_xml(URL)

  x <- as.integer(node2text(myXML, "TaxId"))
  y <- node2text(myXML, "Organism")

  tID <- data.frame()
  if (length(x) > 0 && length(y) > 0) {
    tID <- data.frame(taxID = x, organism = y)
  }
  return(tID)
}



# ==   3.13  UniProtIDmap()  ===================================================
UniProtIDmap <- function (s, mapFrom = "P_REFSEQ_AC", mapTo = "ACC") {
  # Use UniProt ID mapping service to map one or more IDs
  # Parameters:
  #    s        char  A string of white-space separated IDs
  #    mapFrom  char  the database in which the IDs in s are valid.
  #                     Default is RefSeq protein
  #    mapTo    char  the database in which the target IDs are valid.
  #                     Default is UniProtKB
  # Value
  #    A data frame of mapped IDs, with column names From and To, or an
  #    empty data frame if the mapping was unsuccessful. No rows are returned
  #    for IDs that are not mapped.

  # Initialize curl
  httr::set_config(httr::config(http_version = 0))

  URL <- "https://www.uniprot.org/uploadlists/"
  response <- httr::POST(URL,
                         body = list(from = mapFrom,
                                     to = mapTo,
                                     format = "tab",
                                     query = s))

  if (httr::status_code(response) == 200) { # 200: oK
    myMap <- read.delim(file = textConnection(httr::content(response)),
                        sep = "\t")
    colnames(myMap) <- c("From", "To")
  } else {
    myMap <- data.frame()
    warning(paste("No uniProt ID mapping returned:",
                  "server sent status",
                  httr::status_code(response)))
  }

  return(myMap)
}


# ==   3.14  dbProt2JSON()  ====================================================
dbProt2JSON <- function(thisProt) {
  # Extract all protein related data from myDB and return in JSON format.

  thisData <- list()

  # add a protein table
  sel <- which(myDB$protein$name == thisProt)
  thisData$protein <- myDB$protein[sel, ]

  # add a taxonomy table
  sel <- which(myDB$taxonomy$ID == thisData$protein$taxonomyID)
  thisData$taxonomy <- myDB$taxonomy[sel, ]

  # add the entries for this protein from the  annotation table
  sel <- which(myDB$annotation$proteinID == thisData$protein$ID)
  thisData$annotation <- myDB$annotation[sel, ]
  # our .json convention uses pName and fName as keys, not the db-internal IDs
  # add empty columns for pName and fName
  l <- nrow(thisData$annotation)
  thisData$annotation$pName <- character(l)
  thisData$annotation$fName <- character(l)
  # get the appropriate protein and feature names
  for (i in seq_len(l)) {
    pID <- thisData$annotation$proteinID[i]
    sel <- which(myDB$protein$ID == pID)
    thisData$annotation$pName[i] <- myDB$protein$name[sel]    # store pName
    fID <- thisData$annotation$featureID[i]
    sel <- which(myDB$feature$ID == fID)
    thisData$annotation$fName[i] <- myDB$feature$name[sel]    # store fName
  }

  # add the corresponding feature table
  sel <- which(myDB$feature$ID %in% thisData$annotation$featureID)
  thisData$feature <- myDB$feature[sel, ]

  # remove columns that are not going into JSON output
  thisData$protein$ID           <- NULL
  thisData$annotation$ID        <- NULL
  thisData$annotation$proteinID <- NULL
  thisData$annotation$featureID <- NULL
  thisData$feature$ID           <- NULL

  # create JSON-formatted output
  # ( jsonlite::prettify() is too  wordy for a compact Wikipage )

  out <- character()
  out <- c(out, '{')

  out <- c(out, '  "protein": {')
  sel <- colnames(thisData$protein) != "sequence"
  out <- c(out, sprintf("    %s,", dbRow2JSON(thisData$protein[1, sel],
                                              coll = ",\n    ")))
  out <- c(out, dbSeq2JSON(thisData$protein$sequence[1]))
  out <- c(out, '  },')

  out <- c(out, '  "taxonomy": {')
  out <- c(out, sprintf("    %s", dbRow2JSON(thisData$taxonomy)))
  out <- c(out, '  },')

  out <- c(out, '  "annotation": [')
  for (i in seq_len(nrow(thisData$annotation))) {
    out <- c(out, sprintf("    {%s},", dbRow2JSON(thisData$annotation[i, ])))
  }
  out[length(out)] <- gsub(",$", "", out[length(out)]) # remove last ","
  out <- c(out, '  ],')

  out <- c(out, '  "feature": [')
  sel <- colnames(thisData$feature) != "description"
  for (i in seq_len(nrow(thisData$feature))) {
    out <- c(out, sprintf("    {%s,",
                          dbRow2JSON(thisData$feature[i, sel])))
    out <- c(out, sprintf("     %s},",
                          dbRow2JSON(thisData$feature[i, "description",
                                                      drop = FALSE])))
  }
  out[length(out)] <- gsub(",$", "", out[length(out)]) # remove last ","
  out <- c(out, '  ]')

  out <- c(out, '}')

  return(paste0(out, collapse = "\n"))
}


# ==   3.15  dbSeq2JSON()  =====================================================

dbSeq2JSON <- function(s, nIndents = 4, width = 70) {
  # Turn a sequence into a JSON key-value pair, with the value being a JSON
  # array of elements not exceeding a width of "width", and an indent of
  # "indents" spaces.
  ind <- paste0(rep(" ", nIndents), collapse = "")

  out <- character()
  out <- c(out, sprintf("%s\"sequence\": [", ind))

  for (i in seq_along(s)) {
    l <- nchar(s[i])
    if (l <= width) {
      out <- c(out, s[i])
    } else {
      starts <- seq(1, l, by = width)
      ends <- seq(width, l, by = width)
      if (length(ends) < length(starts)) { ends <- c(ends, l) }
      out <- c(out, sprintf("%s  \"%s\",", ind, substring(s[i], starts, ends)))
    }
  }
  out[length(out)] <- gsub(",$", "", out[length(out)]) # remove last ","

  out <- c(out, sprintf("%s]", ind))
  return(paste0(out, collapse = "\n"))
}


# ==   3.16  dbRow2JSON()  =====================================================

dbRow2JSON <- function(df, coll = ", ") {
  # Turn a single dataframe row into JSON key value pairs, where the keys are the
  # column names. Respects character / numeric mode.
  out <- character()
  for (i in 1:ncol(df)) {
    if (class(df[1, i]) == "integer") {
      val <- sprintf("%d", df[1, i])
    } else if (class(df[1, i]) == "numeric") {
      val <- sprintf("%f", df[1, i])
    } else {
      val <- sprintf("\"%s\"", as.character(df[1, i]))
    }
    out <- c(out, sprintf("\"%s\": %s", colnames(df)[i], val))
  }
  return(paste0(out, collapse = coll))
}


# =    4  TESTS  ===============================================================

if (doTESTS) {
  if (! requireNamespace("testthat", quietly = TRUE)) {
    install.packages("testthat")
  }

  # ToDo: test everything here

}


# [END]
