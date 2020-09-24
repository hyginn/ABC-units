# tocID <- "scripts/ABC-dbUtilities.R"
#
# database utilities for ABC learning units
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                             Line
#TOC> -------------------------------------------------
#TOC>   1        PACKAGES                            32
#TOC>   2        FUNCTIONS                           50
#TOC>   2.01       dbSanitizeSequence()              53
#TOC>   2.02       dbConfirmUnique()                 88
#TOC>   2.03       dbInit()                         106
#TOC>   2.04       dbAutoincrement()                147
#TOC>   2.05       dbAddProtein()                   160
#TOC>   2.06       dbAddFeature()                   180
#TOC>   2.07       dbAddTaxonomy()                  199
#TOC>   2.08       dbAddAnnotation()                215
#TOC>   2.09       dbFetchUniProtSeq()              243
#TOC>   2.10       dbFetchPrositeFeatures()         267
#TOC>   2.11       node2text()                      311
#TOC>   2.12       dbFetchNCBItaxData()             323
#TOC>   2.13       UniProtIDmap()                   362
#TOC>   3        TESTS                              401
#TOC> 
#TOC> ==========================================================================


# =    1  PACKAGES  ============================================================


if (! requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}


if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}


if (! requireNamespace("xml2", quietly = TRUE)) {
  install.packages("xml2")
}


# =    2  FUNCTIONS  ===========================================================


# ==   2.01  dbSanitizeSequence()  =============================================
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


# ==   2.02  dbConfirmUnique()  ================================================
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


# ==   2.03  dbInit()  =========================================================
dbInit <- function() {
  # Return an empty instance of the protein database
  # Open the link and study the schema:
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


# ==   2.04  dbAutoincrement()  ================================================
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


# ==   2.05  dbAddProtein()  ===================================================
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
                    sequence    = dbSanitizeSequence(jsonDF$sequence[i]))
    db$protein <- rbind(db$protein, x)
  }
  return(db)
}


# ==   2.06  dbAddFeature()  ===================================================
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
                    accession   = jsonDF$accession[i])
    db$feature <- rbind(db$feature, x)
  }
  return(db)
}


# ==   2.07  dbAddTaxonomy()  ==================================================
dbAddTaxonomy <- function(db, jsonDF) {
  # Add one or more taxonomy entries to the database db.
  # Parameters:
  #     db      list         A database created with dbInit()
  #     jsonDF  data frame   Taxonomy data imported into a data frame with
  #                            fromJSON()
  for (i in seq_len(nrow(jsonDF))) {
    x <- data.frame(
      ID =  jsonDF$ID[i],
      species = jsonDF$species[i])
    db$taxonomy <- rbind(db$taxonomy, x)
  }
  return(db)
}

# ==   2.08  dbAddAnnotation()  ================================================
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
                    end       = as.integer(jsonDF$end[i]))
    db$annotation <- rbind(db$annotation, x)
  }
  return(db)
}


# ==   2.09  dbFetchUniProtSeq()  ==============================================
dbFetchUniProtSeq <- function(ID) {
  # Fetch a protein sequence from UniProt.
  # Parameters:
  #     ID   char   a UniProt ID (accession number)
  # Value:
  #     char        the sequence
  # If the operation is not successful, a 0-length string is returned

  URL <- sprintf("http://www.uniprot.org/uniprot/%s.fasta", ID)

  response <- httr::GET(URL)

  mySeq <- character()
  if (httr::status_code(response) == 200) {
    x <- as.character(response)
    x <- strsplit(x, "\n")
    mySeq <- dbSanitizeSequence(x)
  }

  return(mySeq)
}


# ==   2.10  dbFetchPrositeFeatures()  =========================================
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

    patt <- sprintf("\\|%s\\|", UniProtID)
    lines <- lines[grep(patt, lines)]

    for (line in lines) {
      tokens <- unlist(strsplit(line, "\\t|\\|"))
      myFeatures <- rbind(myFeatures,
                          data.frame(uID   =  tokens[2],
                                     start =  as.numeric(tokens[4]),
                                     end   =  as.numeric(tokens[5]),
                                     psID  =  tokens[6],
                                     psName = tokens[7]))
    }
  }
  return(myFeatures)
}


# ==   2.11  node2text()  ======================================================
node2text <- function(doc, tag) {
  # an extractor function for the contents of elements
  # between given tags in an XML response.
  # Contents of all matching elements is returned in
  # a vector of strings.
  path <- paste0("//", tag)
  nodes <- xml2::xml_find_all(doc, path)
  return(xml2::xml_text(nodes))
}


# ==   2.12  dbFetchNCBItaxData()  =============================================
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



# ==   2.13  UniProtIDmap()  ===================================================
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


# =    3  TESTS  ===============================================================

if (FALSE) {
  if (! requireNamespace("testthat", quietly = TRUE)) {
    install.packages("testthat")
  }

  # ToDo: test everything here

}


# [END]
