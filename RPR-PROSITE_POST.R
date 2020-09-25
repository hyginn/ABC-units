# tocID <- "RPR-PROSITE_POST.R"
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the RPR-Scripting_data_downloads unit.
#
# Version:  1.2
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    2020 Maintenance
#           1.1    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#           1.0.1  Updates for slightly changed interfaces
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
#TOC>   Section  Title                                                 Line
#TOC> ---------------------------------------------------------------------
#TOC>   1        Constructing a POST command from a Web query            43
#TOC>   1.1        Task - fetchPrositeFeatures() function               148
#TOC>   2        Task solutions                                         156
#TOC> 
#TOC> ==========================================================================


# =    1  Constructing a POST command from a Web query  ========================


if (! requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
# Package information:
#  library(help = httr)       # basic information
#  browseVignettes("httr")    # available vignettes
#  data(package = "httr")     # available datasets




# We have reverse engineered the Web form for a ScanProsite request, and can
# construct a valid POST request from knowing the required field names. The POST
# command is similar to GET(), but we need an explicit request body that
# contains a list of key/value pairs

UniProtID <- "P39678"

URL <- "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"

response <- httr::POST(URL,
                       body = list(meta = "opt1",
                                   meta1_protein = "opt1",
                                   seq = UniProtID,
                                   skip = "on",
                                   output = "tabular"))

# Send off this request, and you should have a response in a few
# seconds. Let's check the status first:

httr::status_code(response)  # If this is not 200, something went wrong and it
                             # makes no sense to continue. If this persists, ask
                             # on the Discussion Board what to do.


# The text contents of the response is available with the
# content() function:
httr::content(response, "text")

# ... should show you the same as the page contents that you have seen in the
# browser. Now we need to extract the data from the page. For this simple
# example we can get away with using regular expressions, but in general we need
# a real XML parser to parse HTML. We'll cover that in a later unit. Here, we
# strsplit() the response into individual lines, since each of our data elements
# is on its own line, and then capture the contents. The way Prosite has
# formatted their HTML we can simply split on the "\\n" newline character - but
# they could write the same valid HTML without any newline-characters at all.
# Understand that we are working with a bit of a "hack" here: exploting
# empirical assumptions rather than a formal specification. But sometimes quick
# and dirty is fine, because quick.

lines <- unlist(strsplit(httr::content(response, "text"), "\\n"))
head(lines)

# Now we define a query pattern for the lines we want:
# we can use the uID, bracketed by two "|" pipe
# characters:

patt <- sprintf("\\|%s\\|", UniProtID)

# ... and select only the lines that match this
# pattern:

( lines <- lines[grep(patt, lines)] )

# ... captures the three lines of output.

# Now we break the lines apart into tokens: this is another application of
# strsplit(), but this time we split either on "pipe" characters, "|" OR on tabs
# "\t". Look at the regex "\\t|\\|" in the strsplit() call:

unlist(strsplit(lines[1], "\\t|\\|"))

# Its parts are (\\t)=tab (|)=or (\\|)=pipe. Both "t" and "|" need to be escaped
# with a backslash. "t" has to be escaped because we want to match a tab (\t),
# not the literal character "t". And "|" has to be escaped because we mean the
# literal pipe character, not its metacharacter meaning OR. Thus sometimes the
# backslash turns a special meaning off, and sometimes it turns a special
# meaning on. Unfortunately there's no easy way to tell - you just need to
# remember the characters - or have a reference handy. The metacharacters are
# (){}[]^$?*+.|&-   ... and some of them have different meanings depending on
# where in the regex they are.

# Let's put the tokens into named slots of a data frame

features <- data.frame()
for (line in lines) {
  tokens <- unlist(strsplit(line, "\\t|\\|"))
  features <- rbind(features,
                    data.frame(uID   =  tokens[2],
                               start =  as.numeric(tokens[4]),
                               end   =  as.numeric(tokens[5]),
                               psID  =  tokens[6],
                               psName = tokens[7],
                               psSeq  = tokens[11]))
}
features

#  This forms the base of a function that collects the features automatically
#  from a PrositeScan result. You can write this!


# ==   1.1  Task - fetchPrositeFeatures() function  ============================


# Task: write a function that takes as input a UniProt ID, fetches the
# features it contains from ScanProsite and returns a data frame as given above, or
# an empty data frame if there is an error.


# =    2  Task solutions  ======================================================


# I have placed such a function into the ABC-dbUtilities.R script: look it up by
# clicking on  dbFetchPrositeFeatures() in the Environment pane.

# Test:
dbFetchPrositeFeatures("Q5KMQ9")




# [END]
