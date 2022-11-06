# tocID <- "BIN-Data_integration.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the BIN-Data_integration unit.
#
# Version:  2.0
#
# Date:     2018 - 10  -  2022 - 11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           2.0    2022 Complete rewrite since UniProt Rest API changed
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
#TOC>   Section  Title                                         Line
#TOC> -------------------------------------------------------------
#TOC>   1        Identifier mapping                              48
#TOC>   2        Helper function                                 95
#TOC>   3        Stepping through the process                   145
#TOC>   4        Inspect the results                            222
#TOC>   5        Digression: Cross-referencing tables           249
#TOC>   6        Appendix: Sequence clusters                    298
#TOC> 
#TOC> ==========================================================================


# =    1  Identifier mapping  ==================================================

# Identifier mapping is still a incompletley solved issue of database
# interoperability, as most of you realized when the were asked to get the
# UniProt ID of a protein they found in RefSeq. Why does the NCBI record not
# include this information? Why?

# So you need to spend effort, time, and resources to obtain it yourself.

# UniProt used to provide a well-designed ID mapping tool, but it has mutated
# into ... a thing. What used to be a small, lightweight tool for data
# cross-referencing can now do so very much more. But not as easily. What used
# to be as simple as sending a single string via a Web-browser request is now a
# process that needs several steps. I am sure there are important reasons for
# the massively increased complexity ... but. Let me just say: I'm not a fan. It
# works - but it is neither easy to understand nor instructive to a newcomer.
# And it is not elegant. Too bad. But it works.

# Anyway: here is what needs to happen.
#   - Send a request to the UniProt API
#   - Capture a response with a job ID for your request
#   - Check in from time to time whether the UniProt server completed
#       your job
#   - If yes, get the URL where your results can be found
#   - Retrieve the results
#   - Parse the results into a format you can read.

# The UniProt Web API for the ID mapping service supports a "RESTful API":
# responses can be obtained via a browser request. Such requests are commonly
# sent via the GET or POST verbs that a web server responds to, when a client
# asks for data. GET requests are visible in the URL of the request; POST
# requests are not directly visible, they are commonly used to send the contents
# of forms, or when transmitting larger, complex data items. The UniProt ID
# mapping sevice can accept long lists of IDs, thus using the POST mechanism
# makes sense. GET() and  POST() functions are available as part of the httr::
# package.

# To begin, we make sure httr:: is available. The package supports sending and
# receiving data via the http protocol, just like a Web browser.
if (! requireNamespace("httr", quietly=TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes - useful!


# =    2  Helper function  =====================================================

# The code that appears below originated with UniProt at:
# https://www.uniprot.org/help/id_mapping. I adapted a few parts for style, and
# added comments, and refactored here and there. I also renamed a function and
# most variables (cf. http://thecodelesscode.com/case/220 ). A function was
# called isJobReady(), but it returns TRUE or FALSE. So that's the wrong name:
# it needs to be called jobIsReady(). Right? See, the code should NOT read:

# if (isJobReady()) {  ...
#
#    but it should read:
#
# if (jobIsReady(jobID)) { ...
#
# and actually, since the parameter is already called "jobID", an even better
# name is just isReady() - as in :
#
# if (isReady(jobID)) { ...

# Also: refactoring: a function called isReady() should NOT go and hang out in a
# waiting loop, moreover without feedback to the user! Like, what were you
# thinking? It should do what it's named for because that's what you expect it
# to do. Waiting, and checking back needs to be handled differently. It's
# called: "separation of concerns" and we do that for a reason. So I changed all
# that.


isReady <- function(jobId) {
  # Check whether a given jobID has completed, and if so, whether there
  # were errors or a result.

  URL <- sprintf("https://rest.uniprot.org/idmapping/status/%s", jobId)
  r <- httr::GET(url = URL, httr::accept_json())
  status <- httr::content(r, as = "parsed")

  if (!is.null(status[["messages"]])) {  # Input errors or similar
    stop(status[["messages"]])
  }

  if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
    ready <- TRUE
  } else {
    ready <- FALSE
  }

  return(ready)
}


# =    3  Stepping through the process  ========================================


# We will walk through the process with the refSeqID of yeast Mbp1 and Swi4, and
# we will also enter a dummy ID to check what happens if the ID can't be mapped:
myQueryIDs <- c("NP_010227", "NP_00000", "NP_011036")


# Prepare the request
# ===================
myTask <- list(
  ids = paste(myQueryIDs, collapse =","),
  from = "RefSeq_Protein",
  to = "UniProtKB")


# First interaction:
# ==================
#    send the request,  save the response
r <- httr::POST(url = "https://rest.uniprot.org/idmapping/run/",
                body = myTask,
                encode = "multipart",
                httr::accept_json())

(serverResponse <- httr::content(r, as = "parsed"))

# If everything is sane, this should return a jobID



# Second interaction:
# ===================
#    Check if the results are available

if (isReady(serverResponse$jobId)) {
  cat("Job results available. Please proceed.\n\n")
} else {
  cat("Job results not yet available. Please try again.\n\n")
}



# Third interaction:
# ==================
#    Retrieve the completed job details

URL <- sprintf("https://rest.uniprot.org/idmapping/details/%s",
               serverResponse$jobId)
r <- httr::GET(url = URL, httr::accept_json())
(jobDetails <- httr::content(r, as = "parsed"))



# Fourth interaction:
# ===================
#    Retrieve the actual results at the location specified by the "redirectURL"

#    Where does the URL redirect to? Could be an idmapping stream, or a
#    results stream ...
URL <- jobDetails$redirectURL

if (grepl("/idmapping/results/", URL)) {
  URL <- gsub("/idmapping/results/", "/idmapping/stream/", URL)
} else {
  URL <- gsub("/results/", "/results/stream/", URL)
}

#    We'll request results in a TSV format. Cf.
#    https://www.uniprot.org/help/api_queries#what-formats-are-available
#    Then we read the result into a table.

URL <- paste(URL, "?format=tsv", sep = "")

r <- httr::GET(url = URL, httr::accept_json())
resultsTable <- read.delim(text = httr::content(r))


# =    4  Inspect the results  =================================================

resultsTable[ , 1:3]

# If this works as expected, you should see:
#        From  Entry Entry.Name
# 1 NP_010227 P39678 MBP1_YEAST
# 2 NP_011036 P25302 SWI4_YEAST

# ... and note that there are only two entries, because nothing was returned
# for the dummy "RefSeq ID" NP_00000
#
# Not even an error.

# I have wrapped this into a function - UniProtIDmap() - and added it to the
# ABC-dbUtilities.R script. It is loaded into your workspace during normal
# startup of the course project. Try

?UniProtIDmap

UniProtIDmap(c("NP_010227",
               "NP_011036",
               "NP_012881",
               "NP_013729",
               "NP_012165"))


# =    5  Digression: Cross-referencing tables  ================================

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



# =    6  Appendix: Sequence clusters  =========================================

# This is based on the example that was posted on the UniProt help page: mapping two UniProt IDs to their UniRef90 cluster of similar sequences.

myTask <- list(
  ids = "P21802,P12345",
  from = "UniProtKB_AC-ID",
  to = "UniRef90"
)

URL <- "https://rest.uniprot.org/idmapping/run/"

r <- httr::POST(url = URL,
                body = myTask,
                encode = "multipart",
                httr::accept_json())

serverResponse <- httr::content(r, as = "parsed")

# Wait until this is ready
if (isReady(serverResponse$jobId)) {
  cat("Job results available. Please proceed.\n\n")
} else {
  cat("Job results not yet available. Please try again.\n\n")
}

# Get the results
URL <- sprintf("https://rest.uniprot.org/idmapping/details/%s",
               myResponse$jobId)
r <- httr::GET(url = URL, httr::accept_json())
myDetails <- httr::content(r, as = "parsed")
URL <- myDetails$redirectURL

if (grepl("/idmapping/results/", URL)) {
  URL <- gsub("/idmapping/results/", "/idmapping/stream/", URL)
} else {
  URL <- gsub("/results/", "/results/stream/", URL)
}

URL <- paste(URL, "?format=tsv", sep = "")
r <- httr::GET(url = URL, httr::accept_json())
resultsTable <- read.table(text = httr::content(r),
                           sep = "\t",
                           header=TRUE)



# [END]
