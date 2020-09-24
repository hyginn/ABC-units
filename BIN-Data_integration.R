# tocID <- "BIN-Data_integration.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-Data_integration unit.
#
# Version:  1.2
#
# Date:     2018-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2020 Maintenance and updates
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
#           1.0.1  Bugfix: UniProt ID Mapping service API change
#           1.0    First live version
#
#
# TODO:
#           Develop a fungi-specific BioMart example.
#           (cf.
# https://cran.r-project.org/web/packages/biomartr/vignettes/Functional_Annotation.html )
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                             Line
#TOC> -------------------------------------------------
#TOC>   1        Identifier mapping                  42
#TOC>   2        Cross-referencing tables           165
#TOC>
#TOC> ==========================================================================


# =    1  Identifier mapping  ==================================================

# UniProt provides a well-designed ID mapping tool that can be accessed
# online at     http://www.uniprot.org/mapping/
#
# Here we will use the UniProt Web API for this tool to map identifiers. The
# UniProt ID mapping service supports a "RESTful API": responses can be obtained
# simply via a Web- browsers request. Such requests are commonly sent via the
# GET or POST verbs that a Webserver responds to, when a client asks for data.
# GET requests are visible in the URL of the request; POST requests are not
# directly visible, they are commonly used to send the contents of forms, or
# when transmitting larger, complex data items. The UniProt ID mapping sevice
# can accept long lists of IDs, thus using the POST mechanism makes sense. GET()
# and  POST() functions are part of the httr package.

# To begin, we load  httr, which supports sending and receiving data via the
# http protocol, just like a Web browser.
if (! requireNamespace("httr", quietly=TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes
#  data(package = "httr")     # available datasets


# We will walk through the process with the refSeqID
# of yeast Mbp1 and Swi4, and we will also enter a dummy ID to check what
# happens if the ID can't be mapped:
myQueryIDs <- "NP_010227 NP_00000 NP_011036"


# The UniProt ID mapping service API is very straightforward to use: just define
# the URL of the server and send a list of items labelled as "query" in the body
# of the request. GET() and POST() are functions from httr.

# Note. A recent bug in the interaction between the server expectations and the
# curl client libraries requires the following initialization
httr::set_config(httr::config(http_version = 0))
# cf. https://stackoverflow.com/questions/44610845/stream-error-in-the-http-2-framing-layer-bigrquery-commands-error-in-r-studio-b


URL <- "https://www.uniprot.org/mapping/"
response <- httr::POST(URL,
                       body = list(from = "P_REFSEQ_AC",   # Refseq Protein
                                   to = "ACC",             # UniProt ID
                                   format = "tab",
                                   query = myQueryIDs))

cat(httr::content(response))

# We need to check the status code - if it is not 200, an error ocurred and we
# can't process the result:
httr::status_code(response)

# If the query is successful, tabbed text is returned. We can assign that to a
# data frame. Note that we use textConnection() to read data directly from a char object, which can go in the spot where read.delim() expects a file-name argument.

myMappedIDs <- read.delim(file = textConnection(httr::content(response)),
                          sep = "\t",
                          stringsAsFactors = FALSE)
myMappedIDs

# If this works as expected, you should see:
#        From     To
# 1 NP_010227 P39678
# 2 NP_011036 P25302
#
# ... and note that there are only two entries, because nothing was returned
# for the dummy "RefSeq ID" NP_00000

# If the query can't be fulfilled because of a problem with the server, a
# WebPage is returned. But the server status is also returned and we can check
# the status code. I have lately gotten many "503" status codes: Server Not
# Available...

# We wrap this into a function:

myIDmap <- function (s, mapFrom = "P_REFSEQ_AC", mapTo = "ACC") {
  # Use UniProt ID mapping service to map one or more IDs
  # Parameters:
  #    s  char  A string of separated IDs
  #    mapFrom  char  the database in which the IDs in s are valid. Default
  #                     is RefSeq protein
  #    mapTo    char  the database in which the target IDs are valid. Default
  #                     is UniProtKB
  # Value
  #    a data frame of mapped IDs, with column names From and To, or an
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
                        sep = "\t",
                        stringsAsFactors = FALSE)
    colnames(myMap) <- c("From", "To")
  } else {
    myMap <- data.frame()
    warning(paste("No uniProt ID mapping returned:",
                  "server sent status",
                  httr::status_code(response)))
  }

  return(myMap)
}

# Try it out ...
myIDmap("NP_010227 NP_011036 NP_012881 NP_013729 NP_012165")

# A function UniProtIDmap() is in the ABC-dbUtilities.R script and it is loaded
# into your workspace on startup.


# =    2  Cross-referencing tables  ============================================

# Sometimes we get the IDs we need to map in a large table, e.g. from a list of
# genes in a model organism database such as SGD, or from the Human Genen
# Nomenclature commission. How do we map one set of identifiers to another one?

# The function to use is match().
# Here is a tiny set of identifiers taken from a much larger table to
# illustrate the principle:
#

myIDs <- data.frame(uID =   c("P38903", "P31383", "P47177", "P47096", "Q07747",
                              "Q08641", "P47129", "P52910", "P00330", "P81450"),
                    name =  c("2A5D", "2AAA", "2NDP", "3HAO", "AAD4",
                              "AB140", "ACF4", "ACS2", "ADH1", "ATP18"),
                    refID = c("NP_014657", "NP_009386",
                              "NP_012683", "NP_012559",
                              "NP_010038", "NP_014882",
                              "NP_012616", "NP_013254",
                              "NP_014555", "NP_013629"))

myIDs

# Say we want to map "NP_010038", "NP_012559", and "NP_013629", in that order to
# their gene names.
myQuery <- c("NP_010038", "NP_999999", "NP_013629")

# %in% will only tell us if these IDs are present in the table:
myQuery %in% myIDs$refID

# ... but not where they are located. But match() does what we need here:
match(myQuery, myIDs$refID)

# ... and we can use the result to subset the column that we want to map to:
myIDs$name[match(myQuery, myIDs$refID)]

# Note that the output preserves the NA - i.e. the length of the mapped
# values is exactly the same as the length of the query.

# task: map the three genes to their UniProt Identifier.


#
# Note: if you want to do very many queries in very large tables, use the
# fmatch() function in the "fastmatch" package for a considerable
# speedup.




# [END]
