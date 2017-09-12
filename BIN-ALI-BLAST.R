# BIN-ALI-BLAST.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-ALI-BLAST unit.
#
# Version:  0.1
#
# Date:     2017  08  28
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.1    First code copied from 2016 material.

#
# TODO:
#
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...

# ==============================================================================

# = 1 ___Section___

# BLAST.R
#
# Purpose: Send off one BLAST search and return parsed list of results
#          This script uses the BLAST URL-API
#          (Application Programming Interface) at the NCBI.
#          Read about the constraints here:
# http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYP=DeveloperInfo
#
#
# Version: 1.0
# Date:    2016-09
# Author:  Boris Steipe
#
#
# ToDo:
# Notes:   The bioconducter "annotate" package contains code for BLAST searches,
#          in case you need to do something more involved.
#
# ==============================================================================


# Dependencies:  myEmail must exist as a global variable with
#                     your valid email adress
#                waitTimer() must be loaded (it should have been loaded from
#                     .utilities.R, which was sourced via .Rprofile)


# library to interface with WebServers and process their XML/HTML
# responses
if (!require(xml2, quietly = TRUE)) {
  install.packages("xml2")
  library(xml2)
}

if (!require(httr, quietly = TRUE)) {
  install.packages("httr")
  library(httr)
}



parseBLAST_XML <- function(hit) {
  # parse one BLAST hit XML node with the xml2 package;
  # return a list

  h <- list()
  h$id <-  xml_text(xml_find_first(hit, ".//Hit_accession"))
  h$def <- xml_text(xml_find_first(hit, ".//Hit_def"))
  h$bestE <- Inf
  h$sumLen <- 0
  h$sumId <- 0
  h$sumGap <- 0
  hsps <- xml_find_all(hit, ".//Hsp")
  h$Hsp <- list()
  h$nHsps <- length(hsps)
  if (h$nHsps > 0) {
    for (i in 1:length(hsps)) {
      h$Hsp[[i]] <- list()
      h$Hsp[[i]]$e <-          xml_numeric(hsps[i], ".//Hsp_evalue")
      h$Hsp[[i]]$q_from <-     xml_numeric(hsps[i], ".//Hsp_query-from")
      h$Hsp[[i]]$q_to <-       xml_numeric(hsps[i], ".//Hsp_query-to")
      h$Hsp[[i]]$h_from <-     xml_numeric(hsps[i], ".//Hsp_hit-from")
      h$Hsp[[i]]$h_to <-       xml_numeric(hsps[i], ".//Hsp_hit-to")
      h$Hsp[[i]]$h_identity <- xml_numeric(hsps[i], ".//Hsp_identity")
      h$Hsp[[i]]$h_gaps <-     xml_numeric(hsps[i], ".//Hsp_gaps")
      h$Hsp[[i]]$h_len <-      xml_numeric(hsps[i], ".//Hsp_align-len")
      h$Hsp[[i]]$qseq <- xml_text(xml_find_first(hsps[i], ".//Hsp_qseq"))
      h$Hsp[[i]]$mid <-  xml_text(xml_find_first(hsps[i], ".//Hsp_midline"))
      h$Hsp[[i]]$hseq <- xml_text(xml_find_first(hsps[i], ".//Hsp_hseq"))

      h$bestE <- min(h$bestE, h$Hsp[[i]]$e)
      h$sumLen <- h$sumLen + h$Hsp[[i]]$h_len
      h$sumId <- h$sumId + h$Hsp[[i]]$h_identity
      h$sumGap <- h$sumGap + h$Hsp[[i]]$h_gaps
    }
  }
  return(h)
}

xml_numeric <- function(n, p) {
  # Utility: return first node matching xpath p in XML node n as numeric
  return(as.numeric(xml_text(xml_find_first(n, p))))
}


BLAST <- function(q,
                  db = "refseq_protein",
                  nHits = 30,
                  E = 3,
                  limits = "\"\"",
                  email = myEMail,
                  rid = "",
                  quietly = FALSE) {
  # Purpose:
  #     Basic BLAST search
  # Version: 1.0
  # Date:    2016-09
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
  #     email: a valid email address, defaults to global value myEMail
  #     quietly: controls printing of wait-time progress bar
  # Value:
  #     result: list of resulting hits and some metadata

  results <- list()
  results$rid <- rid
  results$rtoe <- 0

  if (rid == "") {  # prepare, send and analyse query
    results$query <- paste(
      "https://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
      "?",
      "QUERY=", q,
      "&DATABASE=", db,
      "&HITLIST_SIZE=", as.character(nHits),
      "&EXPECT=", as.character(E),
      "&PROGRAM=", "blastp",
      "&ENTREZ_QUERY=", limits,
      "&NOHEADER=", "true",
      "&EMAIL=", email,
      "&CMD=Put",
      sep = "")

    # send it off ...
    response <- read_xml(results$query, as_html = TRUE)

    # find the comment node that contains the information we need
    # using an xpath expression
    info <- xml_find_first(response,
                           "//comment()[contains(., \"QBlastInfo\")]")

    info <- xml_text(info)  # extract its contents

    # parse
    results$rid  <- regmatches(info,
                               regexec("RID = (\\w+)",  info))[[1]][2]
    results$rtoe <- regmatches(info,
                               regexec("RTOE = (\\d+)", info))[[1]][2]
    results$rtoe <- as.numeric(results$rtoe)
  } # done analysing query

  # Now we wait ...
  if (quietly) {
    Sys.sleep(results$rtoe)
  } else {
    cat(sprintf("BLAST is processing %s:\n", results$rid))
    waitTimer(results$rtoe)
  }

  # retrieve results from BLAST server
  URL <- paste("https://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
               "?",
               "RID=", results$rid,
               "&FORMAT_TYPE=", "XML",
               "&EMAIL=", email,
               "&CMD=Get",
               sep = "")
  raw <- GET(URL)

  timeOut <- 300
  nWait <- 0
  while (raw$headers["content-type"] == "text/html" && nWait <= (timeOut/10)) {
    cat("Doesn't seem to be done. Wait some more (or click STOP to abort)\n")
    waitTimer(10)
    nWait <- nWait + 1
    raw <- GET(URL)
  }

  # If we get to here, we received some result. But what?
  if (raw$headers["content-type"] == "text/html") { # Still HTML? Didn't complete ...
    stop(sprintf("Query >>%s<< didn't complete.", results$rid))
  } else if (raw$headers["content-type"] == "application/xml") { # Good!
    response <- read_xml(raw)
  } else { # Unknown, abort.
    stop(sprintf("Unknown response type: >>%s<<.", raw$headers["content-type"]))
  }

  hits <- xml_find_all(response, ".//Hit")

  if (length(hits) == 0) {
    s <- "No hit returned.\n"
    s <- c(s, sprintf("Check your query string:\n>>%s<<\n", results$query))
    s <- c(s, sprintf("and/or try again later by typing:\n", results$rid))
    s <- c(s, sprintf("   BLAST(rid = \"%s\")\n", results$rid))
    stop(paste(s, collapse = ""))
  }

  results$hits <- list()

  for (i in 1:length(hits)) {
    results$hits[[i]] <- parseBLAST_XML(hits[i])
  }

  return(results)
}



# = 1 Tasks




# ==== TESTS ===================================================================

# q   <- paste("IYSARYSGVDVYEFIHSTGSIMKRKKDDWVNATHI",   # Mbp1 APSES domain
#              "LKAANFAKAKRTRILEKEVLKETHEKVQGGFGKYQ",
#              "GTWVPLNIAKQLAEKFSVYDQLKPLFDFTQTDGSASP",
#              sep="")
# q <- "NP_010227"
# fungi <- "txid4751[ORGN]"
#
# test <- BLAST("NP_010227",
#               nHits = 1000,
#               E = 0.01,
#               limits = fungi)
# length(test$hits)



# [END]
