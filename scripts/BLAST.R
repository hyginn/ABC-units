# BLAST.R
#
# Purpose: Send off one BLAST search and return parsed list of results
#          This script uses the BLAST URL-API
#          (Application Programming Interface) at the NCBI.
#          Read about the constraints here:
#          https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
#
#
# Version: 3.2
# Date:    2016 09 - 2020 09
# Author:  Boris Steipe
#
# Versions:
#    3.2   2020 updates
#    3.1   Change from require() to requireNamespace(),
#          use <package>::<function>() idiom throughout
#    3.0   parsing logic had not been fully implemented; Fixed.
#    2.1   bugfix in BLAST(), bug was blanking non-split deflines;
#          refactored parseBLASTalignment() to handle lists with multiple hits.
#    2.0   Completely rewritten because the interface completely changed.
#          Code adpated in part from NCBI Perl sample code:
#          $Id: web_blast.pl,v 1.10 2016/07/13 14:32:50 merezhuk Exp $
#    1.0   first version posted for BCH441 2016, based on BLAST - API
#
# ToDo:    Return the organism/strain name in the output, and propagate
#          into MYSPE selection script.
#
# Notes:   This is somewhat pedestrian, but apparently there are currently
#          no R packages that contain such code.
#
# ==============================================================================


if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}


BLAST <- function(Q,
                  db = "refseq_protein",
                  nHits = 30,
                  E = 0.1,
                  limits = "",
                  rid = "",
                  query = "",
                  quietly = FALSE,
                  myTimeout = 120) {
    # Purpose:
    #     Basic BLAST search
    #
    # Parameters:
    #     Q: query - either a valid ID or a sequence
    #     db: "refseq_protein" by default,
    #         other legal values include: "nr", "pdb", "swissprot" ...
    #     nHits: number of hits to maximally return
    #     E: E-value cutoff. Do not return hits whose score would be expected
    #        to occur E or more times in a database of random sequence.
    #     limits: a valid ENTREZ filter
    #     rid: a request ID - to retrieve earlier search results
    #     query: the actual query string (needed when retrieving results
    #            with an rid)
    #     quietly: controls printing of wait-time progress bar
    #     timeout: how much longer _after_ rtoe to wait for a result
    #              before giving up (seconds)
    # Value:
    #     result: list of process status or resulting hits, and some metadata


    EXTRAWAIT <- 10 # duration of extra wait cycles if BLAST search is not done

    results <- list()
    results$query = query
    results$rid <- rid
    results$rtoe <- 0

    if (rid == "") {  # If no rid is available, spawn a search.
                      # Else, proceed directly to retrieval.

      # prepare query, GET(), and parse rid and rtoe from BLAST server response
      results$query <- paste0("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi",
                              "?",
                              "CMD=Put",
                              "&PROGRAM=", "blastp",
                              "&QUERY=", URLencode(Q),
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
      response <- httr::GET(results$query)
      if (httr::http_status(response)$category != "Success" ) {
        stop(sprintf("PANIC: Can't send query. BLAST server status error: %s",
                     httr::http_status(response)$message))
      }

      txt <- httr::content(response, "text", encoding = "UTF-8")

      patt <- "RID = (\\w+)" # match the request id
      results$rid  <- regmatches(txt, regexec(patt,  txt))[[1]][2]

      patt <- "RTOE = (\\d+)" # match the expected completion time
      results$rtoe <- as.numeric(regmatches(txt, regexec(patt, txt))[[1]][2])

      # Now we wait ...
      if (quietly) {
        Sys.sleep(results$rtoe)
      } else {
        cat(sprintf("BLAST is processing %s (rtoe is %d seconds):\n",
                    results$rid,
                    results$rtoe))
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
      response <- httr::GET(checkStatus)
      if (httr::http_status(response)$category != "Success" ) {
        stop(sprintf("PANIC: Can't check status. BLAST server status error: %s",
                     httr::http_status(response)$message))
      }

      txt <- httr::content(response, "text", encoding = "UTF-8")

      if (length(grep("Status=WAITING",  txt)) > 0) {
        myTimeout <- myTimeout - EXTRAWAIT

        if (myTimeout <= 0) { # abort
          cat("BLAST search not concluded before timeout. Aborting.\n")
          cat(sprintf("%s  BLASThits <- BLAST(rid=\"%s\")\n\n",
                      "Trying checking back later with >\n",
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

    response <- httr::GET(retrieve)
    if (httr::http_status(response)$category != "Success" ) {
      stop(sprintf("PANIC: Can't retrieve. BLAST server status error: %s",
                   httr::http_status(response)$message))
    }

    txt <- httr::content(response, "text", encoding = "UTF-8")

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
    if (length(x) > 0) {
      txt[x-1] <- paste0(txt[x-1], txt[x])
      txt <- txt[-x]
    }

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
  # Parse data from a character vector containing a BLAST hit
  # Parameters:
  #    hit  char   one BLAST hit as char vector
  # Value:
  #          list   $def          chr   defline
  #                 $accession    chr   accession number
  #                 $organism     chr   complete organism definition
  #                 $species      chr   binomial species
  #                 $E            num   E value
  #                 $lengthAli    num   length of the alignment
  #                 $nIdentitites num   number of identities
  #                 $nGaps        num   number of gaps
  #                 $Qbounds      num   2-element vector of query start-end
  #                 $Sbounds      num   2-element vector of subject start-end
  #                 $Qseq         chr   query sequence
  #                 $midSeq       chr   midline string
  #                 $Sseq         chr   subject sequence

  getToken <- function(patt, v) {
    # get the first token identified by pattern patt in character vector v
    v <- v[grep(patt, v)]
    if (length(v) > 1) { v <- v[1] }
    if (length(v) == 0) { token <- NA
    } else {
      token <- regmatches(v, regexec(patt, v))[[1]][2] }
    return(token)
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
  if (length(x) >= 2) {
    h$species <- paste(x[1], x[2])
  } else if (length(x) == 1) {
    h$species <- paste(x[1], "sp.")
  } else {
    h$species <- NA
  }

  # E-value
  h$E <- as.numeric(getToken("Expect\\s*=(.+?), Method", hit))

  # length of alignment
  h$lengthAli <- as.numeric(getToken("^\\s*Length\\s*=(.+)$", hit))

  # number of identities
  h$nIdentities <- as.numeric(getToken("^\\s*Identities\\s*=(.+?)/", hit))

  # number of gaps
  h$nGaps <- as.numeric(getToken("\\s*Gaps\\s*=(.+?)/", hit))

  # split up alignment section
  idx <- grep("^Query ", hit)
  Que <- hit[idx]
  Mid <- hit[idx + 1]
  Sbj <- hit[idx + 2]

  # first and last positions
  h$Qbounds <- c(start = 0, end = 0)
  h$Qbounds[1] <- as.numeric(getToken("^Query\\s*(\\d+)", Que[1]))
  h$Qbounds[2] <- as.numeric(getToken("\\s*(\\d+)\\s*$", Que[length(Que)]))

  h$Sbounds <- c(start = 0, end = 0)
  h$Sbounds[1] <- as.numeric(getToken("^Sbjct\\s*(\\d+)", Sbj[1]))
  h$Sbounds[2] <- as.numeric(getToken("\\s*(\\d+)\\s*$", Sbj[length(Sbj)]))

  # aligned sequences
  for (i in seq_along(Que)) {
    patt <- ("^\\s*Query\\s*\\d+\\s*([A-Za-z-]+)") # capture aligned string
    m <- regexec(patt, Que[i])
    iFirst <- m[[1]][2]
    iLast <- iFirst + attr(m[[1]], which = "match.length")[2] - 1
    Que[i] <- substring(Que[i], iFirst, iLast)
    Mid[i] <- substring(Mid[i], iFirst, iLast)
    Sbj[i] <- substring(Sbj[i], iFirst, iLast)
  }

  h$Qseq   <- paste0(Que, collapse = "")
  h$midSeq <- paste0(Mid, collapse = "")
  h$Sseq   <- paste0(Sbj, collapse = "")

  return(h)
}


# ==== TESTS ===================================================================

if (FALSE) {
  # define query:
  q   <- paste("IYSARYSGVDVYEFIHSTGSIMKRKKDDWVNATHI", # Mbp1 APSES domain
               "LKAANFAKAKRTRILEKEVLKETHEKVQGGFGKYQ",
               "GTWVPLNIAKQLAEKFSVYDQLKPLFDFTQTDGSASP",
               sep="")
  # or ...
  q <- "NP_010227" # refseq ID

  test <- BLAST(q,
                nHits = 100,
                E = 0.001,
                rid = "",
                limits = "txid4751[ORGN]")  # Fungi
  str(test)
  length(test$hits)
}

# [END]

