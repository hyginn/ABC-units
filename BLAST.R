# BLAST.R
#
# Purpose: Send off one BLAST search and return parsed list of results
#          This script uses the BLAST URL-API
#          (Application Programming Interface) at the NCBI.
#          Read about the constraints here:
# https://ncbi.github.io/blast-cloud/dev/api.html
#
#
# Version: 2.0
# Date:    2016 09 - 2017 09
# Author:  Boris Steipe
#
# Versions:
#    2.0   Completely rewritten because the interface completely changed.
#          Code adpated in part from NCBI Perl sample code:
#          $Id: web_blast.pl,v 1.10 2016/07/13 14:32:50 merezhuk Exp $
#
#    1.0   first version posted for BCH441 2016, based on BLAST - API
#
# ToDo:
#
# Notes:   This is somewhat pedestrian, but apparently there are currently
#          no R packages that contain such code.
#
# ==============================================================================


if (!require(httr, quietly = TRUE)) {
  install.packages("httr")
  library(httr)
}


BLAST <- function(q,
                  db = "refseq_protein",
                  nHits = 30,
                  E = 0.1,
                  limits = "",
                  rid = "",
                  quietly = FALSE,
                  myTimeout = 120) {
    # Purpose:
    #     Basic BLAST search
    # Version: 2.0
    # Date:    2017-09
    # Author:  Boris Steipe
    #
    # Parameters:
    #     q: query - either a valid ID or a sequence
    #     db: "refseq_protein" by default,
    #         other legal valuses include: "nr", "pdb", "swissprot" ...
    #     nHits: number of hits to maximally return
    #     E: E-value cutoff. Do not return hits whose score would be expected
    #        to occur E or more times in a database of random sequence.
    #     limits: a valid ENTREZ filter
    #     rid: a request ID - to retrieve earleir search results
    #     quietly: controls printing of wait-time progress bar
    #     timeout: how much longer _after_ rtoe to wait for a result
    #              before giving up (seconds)
    # Value:
    #     result: list of resulting hits and some metadata


    EXTRAWAIT <- 10 # duration of extra wait cycles if BLAST search is not done

    results <- list()
    results$rid <- rid
    results$rtoe <- 0

    if (rid == "") {  # we skip, and proceed directly to retrieval
                      # if rid is not the empty string

      # prepare query, GET(), and parse rid and rtoe from BLAST server response
      results$query <- paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                              "?",
                              "CMD=Put",
                              "&PROGRAM=", "blastp",
                              "&QUERY=", URLencode(q),
                              "&DATABASE=", db,
                              "&MATRIX=", "BLOSUM62",
                              "&EXPECT=", as.character(E),
                              "&HITLIST_SIZE=", as.character(nHits),
                              "&ALIGNMENTS=", as.character(nHits),
                              "&FORMAT_TYPE=Text")

      if (limits != "") {
        results$query <- paste0(
          results$query,
          "&ENTREZ_QUERY=", limits)
      }

      # send it off ...
      response <- GET(results$query)
      if (http_status(response)$category != "Success" ) {
        stop(sprintf("PANIC: Can't send query. BLAST server status error: %s",
                     http_status(response)$message))
      }

      txt <- content(response, "text", encoding = "UTF-8")

      patt <- "RID = (\\w+)" # match the request id
      results$rid  <- regmatches(txt, regexec(patt,  txt))[[1]][2]

      patt <- "RTOE = (\\d+)" # match the expected completion time
      results$rtoe <- as.numeric(regmatches(txt, regexec(patt, txt))[[1]][2])

      # Now we wait ...
      if (quietly) {
        Sys.sleep(results$rtoe)
      } else {
        cat(sprintf("BLAST is processing %s:\n", results$rid))
        waitTimer(results$rtoe)
      }

    } # done sending query and retrieving rid, rtoe

    # Enter an infinite loop to check for result availability
    checkStatus <- paste("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                         "?",
                         "CMD=Get",
                         "&RID=", results$rid,
                         "&FORMAT_TYPE=Text",
                         "&FORMAT_OBJECT=SearchInfo",
                         sep = "")

    while (TRUE) {
      # Check whether the result is ready
      response <- GET(checkStatus)
      if (http_status(response)$category != "Success" ) {
        stop(sprintf("PANIC: Can't check status. BLAST server status error: %s",
                     http_status(response)$message))
      }

      txt <- content(response, "text", encoding = "UTF-8")

      if (length(grep("Status=WAITING",  txt)) > 0) {
        myTimeout <- myTimeout - EXTRAWAIT

        if (myTimeout <= 0) { # abort
          cat("BLAST search not concluded before timeout. Aborting.\n")
          cat(sprintf("You could check back later with rid \"%s\"\n",
                      results$rid))
          return(results)
        }

        if (quietly) {
          Sys.sleep(EXTRAWAIT)
        } else {
          cat(sprintf("Status: Waiting. Wait %d more seconds (max. %d more)",
                      EXTRAWAIT,
                      myTimeout))
          waitTimer(EXTRAWAIT)
          next
        }

      } else if (length(grep("Status=FAILED",  txt)) > 0) {
          cat("BLAST search returned status \"FAILED\". Aborting.\n")
          return(results)

      } else if (length(grep("Status=UNKNOWN",  txt)) > 0) {
          cat("BLAST search returned status \"UNKNOWN\".\n")
          cat("This probably means the rid has expired. Aborting.\n")
          return(results)

      } else if (length(grep("Status=READY",  txt)) > 0) {  # Done

          if (length(grep("ThereAreHits=yes",  txt)) == 0) {  # No hits
            cat("BLAST search ready but no hits found. Aborting.\n")
            return(results)

          } else {
            break  # done ... retrieve search result
          }
      }
    } # end result-check loop

    # retrieve results from BLAST server
    retrieve <- paste("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                      "?",
                      "&CMD=Get",
                      "&RID=", results$rid,
                      "&FORMAT_TYPE=Text",
                      sep = "")

    response <- GET(retrieve)
    if (http_status(response)$category != "Success" ) {
      stop(sprintf("PANIC: Can't retrieve. BLAST server status error: %s",
                   http_status(response)$message))
    }

    txt <- content(response, "text", encoding = "UTF-8")

    # txt contains the whole set of results. Process:

    # First, we strsplit() on linebreaks:
    txt <- unlist(strsplit(txt, "\n"))

    # The alignments range from the first line that begins with ">" ...
    iFirst <- grep("^>", txt)[1]

    # ... to the last line that begins with "Sbjct"
    x <- grep("^Sbjct", txt)
    iLast <- x[length(x)]

    # Get the alignments block
    txt <- txt[iFirst:iLast]

    # Drop empty lines
    txt <- txt[!(nchar(txt) == 0)]

    # A line that ends "]" but does not begin ">" seems to be a split
    # defline ... eg.
    #  [1] ">XP_013349208.1 AUEXF2481DRAFT_695809 [Aureobasidium subglaciale "
    #  [2] "EXF-2481]"
    #  Merge these lines to the preceding lines and delete them.
    #
    x <- which(grepl("]$", txt) & !(grepl("^>", txt)))
    txt[x-1] <- paste0(txt[x-1], txt[x])
    txt <- txt[-x]

    # Special case: there may be multiple deflines when the BLAST hit is to
    # redundant, identical sequences. Keep only the first instance.
    iKeep <- ! grepl("^>", txt)
    x <- rle(iKeep)
    x$positions <- cumsum(x$lengths)
    i <- which(x$lengths > 1 & x$values == FALSE)
    if (length(i) > 0) {
      firsts <- x$positions[i] - x$lengths[i] + 1
      iKeep[firsts] <- TRUE
      txt <- txt[iKeep]
    }

    # After this preprocessing the following should be true:
    # - Every alignment block begins with a defline in which the
    #   first character is ">"
    # - There is only one defline in each block.
    # - Lines are not split.

    # Make a dataframe of first and last indices of alignment blocks
    x <- grep("^>", txt)
    blocks <- data.frame(iFirst = x,
                         iLast  = c((x[-1] - 1), length(txt)))

    # Build the hits list by parsing the blocks
    results$hits <- list()

    for (i in seq_len(nrow(blocks))) {
      thisBlock <- txt[blocks$iFirst[i]:blocks$iLast[i]]
      results$hits[[i]] <- parseBLASTalignment(thisBlock)
    }

    return(results)
}

parseBLASTalignment <- function(hit) {
  # parse one BLAST hit;
  # return a list

  if (length(grep("Length", hit)) > 1) {
    stop("Parsing function can't handle multiple HSPs (yet).")
  }

  h <- list()

  # FASTA defline
  h$def <- hit[1]

  # accesion number (ID), use the first if there are several, separated by "|"
  patt <- "^>(.+?)(\\s|\\|)" # from ">" to space or "|"
  h$accession <-  regmatches(h$def, regexec(patt, h$def))[[1]][2]

  # organism
  patt <- "\\[(.+)]"
  h$organism <-  regmatches(h$def, regexec(patt, h$def))[[1]][2]

  # species
  x <- unlist(strsplit(h$organism, "\\s+"))
  if (length(x) < 2) {
    h$species <- NA
  } else {
    h$species <- paste(x[1:2], collapse = " ")
  }

  # E-value
  x <- hit[grep("Expect\\s*=", hit)]
  patt <- "Expect\\s*=\\s*([0-9.eE\\-]+)" #
  h$E <-  as.numeric(regmatches(x, regexec(patt, x))[[1]][2])

  # length of hit and # identities
  x <- hit[grep("Identities\\s*=", hit)]
  patt <- "Identities\\s*=\\s*([0-9]+)/([0-9]+)"
  m <- regexec(patt, x)
  h$lengthAli   <- as.numeric(regmatches(x, m)[[1]][2])
  h$nIdentities <- as.numeric(regmatches(x, m)[[1]][3])

  # number of gaps
  x <- hit[grep("Gaps\\s*=", hit)]
  patt <- "Gaps\\s*=\\s*([0-9]+)"
  h$nGaps <- as.numeric(regmatches(x, regexec(patt, x))[[1]][2])

  # first and last positions
  iAli <- grep("^Query\\s+", hit)
  h$Qbounds <- getAliBounds(hit[iAli])
  h$Sbounds <- getAliBounds(hit[iAli + 2])

  # aligned sequences

  h$Qseq   <- character()
  h$midSeq <- character()
  h$Sseq   <- character()

  for (i in iAli) {
    patt <- "^Query\\s+[0-9]+\\s*"
    first <- attr(regexec(patt, hit[i])[[1]], "match.length") + 1

    patt <- "\\s*[0-9]*\\s*$"
    last <- regexec(patt, hit[i])[[1]][1] - 1

    h$Qseq   <- paste0(h$Qseq,   substr(hit[i],     first, last))
    h$midSeq <- paste0(h$midSeq, substr(hit[i + 1], first, last))
    h$Sseq   <- paste0(h$Sseq,   substr(hit[i + 2], first, last))
  }

  return(h)
}


getAliBounds <- function(s) {
  # get first and last position from a vector of BLAST alignments s
  # value: numeric vector of first and last position
  patt <- "^(Query|Sbjct)\\s+([0-9]+)\\s"
  first <- as.numeric(regmatches(s[1], regexec(patt, s[1]))[[1]][3])

  patt <- "\\s*([0-9]+)\\s*$"
  last <- as.numeric(regmatches(s[length(s)],
                                regexec(patt, s[length(s)]))[[1]][2])
  return(c (first, last))
}



# ==== TESTS ===================================================================

# define query:
# q   <- paste("IYSARYSGVDVYEFIHSTGSIMKRKKDDWVNATHI", # Mbp1 APSES domain sequence
#              "LKAANFAKAKRTRILEKEVLKETHEKVQGGFGKYQ",
#              "GTWVPLNIAKQLAEKFSVYDQLKPLFDFTQTDGSASP",
#              sep="")
# or ...
# q <- "NP_010227" # refseq ID
#
# test <- BLAST(q,
#               nHits = 100,
#               E = 0.001,
#               rid = "",
#               limits = "txid4751[ORGN]")
# length(test$hits)

# [END]

