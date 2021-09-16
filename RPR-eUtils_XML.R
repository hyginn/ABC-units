# tocID <- "RPR-eUtils_XML.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Scripting_data_downloads unit.
#
# Version:  1.2.1
#
# Date:     2017-10  -  2021-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2.1  2021 Maintenance
#           1.2    2020 Updates
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
#           1.0    First ABC units version
#           0.1    First code copied from 2016 material.
#
#
# TODO:
#
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
#TOC>   Section  Title                                       Line
#TOC> -----------------------------------------------------------
#TOC>   1        Working with NCBI eUtils                      43
#TOC>   1.1        Task - fetchNCBItaxData() function         145
#TOC>   2        Task solutions                               152
#TOC> 
#TOC> ==========================================================================


# =    1  Working with NCBI eUtils  ============================================


# To begin, we load the xml2 package that contains functions
# we need to receive and parse html data. NCBI's eUtils send information in
# XML format so we need to be able to parse XML.
if (! requireNamespace("xml2", quietly=TRUE)) {
  install.packages("xml2")
}
# Package information:
#  library(help = xml2)       # basic information
#  browseVignettes("xml2")    # available vignettes
#  data(package = "xml2")     # available datasets



# We will walk through the process with the refSeqID
# of yeast Mbp1
refSeqID <- "NP_010227"


# First we build a query URL...
eUtilsBase <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


# Then we assemble an URL that will search for get the
# unique, NCBI internal identifier,
# for our refSeqID...
URL <- paste(eUtilsBase,
             "esearch.fcgi?",     # ...using the esearch program
                                  # that finds an entry in an
                                  # NCBI database
             "db=protein",
             "&term=", refSeqID,
             sep="")
# Copy the URL and paste it into your browser to see
# what the response should look like.
URL

# To fetch a response in R, we use the function read_xml()
# with our URL as its argument.
( myXML <- xml2::read_xml(URL) )

# This is XML. We can take the response apart into
# its individual components with the as_list() function.

xml2::as_list(myXML)

# Note how the XML "tree" is represented as a list of
# lists of lists ...
# If we know exactly what element we are looking for,
# we can extract it from this structure:
xml2::as_list(myXML)[["eSearchResult"]][["IdList"]][["Id"]][[1]]

# But this is not very robust, it would break with the
# slightest change that the NCBI makes to their data format -
# and the NCBI changes things A LOT!

# Somewhat more robust is to specify the type of element
# we want - its the text contained in an <Id>...</Id>
# element, and use the XPath XML parsing language to
# retrieve it.

xml2::xml_find_all(myXML, "//Id") # returns a "node set"

xml2::xml_text(xml2::xml_find_all(myXML, "//Id")) # returns the contents
                                                  # of the node set

# We will need to do this more than once, so we write a function
# for it...
node2text <- function(doc, tag) {
  # an extractor function for the contents of elements
  # between given tags in an XML response.
  # Contents of all matching elements is returned in
  # a vector of strings.
  path <- paste0("//", tag)
  nodes <- xml2::xml_find_all(doc, path)
  return(xml2::xml_text(nodes))
}

# using node2text() ...
(GID <- node2text(myXML, "Id"))

# The GI is the pivot for data requests at the
# NCBI.

# Let's first get the associated data for this GI
URL <- paste0(eUtilsBase,
              "esummary.fcgi?",
              "db=protein",
              "&id=",
              GID,
              "&version=2.0")
(myXML <- xml2::read_xml(URL))

(taxID <- node2text(myXML, "TaxId"))
(organism <- node2text(myXML, "Organism"))

#  This forms the base of a function that gets taxonomy data
#  from an Entrez result. You can write this!


# ==   1.1  Task - fetchNCBItaxData() function  ================================

# Task: write a function that takes as input a RefSeq ID, fetches the taxonomy
# information, returns a list with taxID and organism, if the operation is
# successful, or a list of length 0 if there is an error.


# =    2  Task solutions  ======================================================

# I have placed such a function into the dbUtilities script: look it up by
# clicking on  dbFetchNCBItaxData() in the Environment pane.

# Test:
dbFetchNCBItaxData("XP_001837394")

# Expected outout:
# ----------------
# taxID                         organism
# 1 240176 Coprinopsis cinerea okayama7#130


# [END]
