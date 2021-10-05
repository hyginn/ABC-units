# tocID <- "RPR-Pipe.R"
#
# Purpose:  A Bioinformatics Course:
#              Discussing pipe operators.
#
# Version:  1.0
#
# Date:     2021  10
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.0    New code
#
#
# TODO:
#   - find more interesting examples
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
#TOC>   Section  Title                            Line
#TOC> ------------------------------------------------
#TOC>   1        Pipe  Concept                      41
#TOC>   2        Nested Expression                  73
#TOC>   3        magrittr:: Pipe                    78
#TOC>   4        Base R Pipe                        93
#TOC>   5        Intermediate Assignment           108
#TOC>   6        Postscript                        127
#TOC>
#TOC> ==========================================================================


# =    1  Pipe  Concept  =======================================================

# Pipes are actually an awesome idea for any code that implements a workflow -
# a sequence of operations, each of which transforms data in a specialized way.
#
# This principle is familiar from maths: chained functions. If have a function
# y = f(x) and want to use those results as in z = g(y), I can just write
# z = g(f(x))
#
# On the unix command line, pipes were used from the very beginning, implemented
# with the "|" pipe character.
#
# In R, the magrittr package provided the %>% operator, and recently the |>
# operator has been introduced into base R.
#
# However there are alternatives: intermediate assignment, and nested functions
# that have always existed in base R anyway.
#
# Let us look at an example. In writing this, I found out that virtually
# ALL non-trivial examples I came up with don't translate well into this idiom
# at all. It is actually quite limited to simple filtering operations on
# data. A more interesting example might be added in the future, let me know if
# you have a good idea.
#
# A somewhat contrived example is to sort a list of files by the
# length of the file names:

myFiles <- list.files(pattern = "\\.R$")

# nchar() gives the number of characters in a string, order() produces indices
# that map an array to its sorted form.
#
# =    2  Nested Expression  ===================================================

myFiles[order(nchar(myFiles))]


# =    3  magrittr:: Pipe  =====================================================

if (! requireNamespace("magrittr", quietly = TRUE)) {
  install.packages("magrittr")
}
# Package information:
#  library(help = magrittr)       # basic information
#  browseVignettes("magrittr")    # available vignettes
#  data(package = "magrittr")     # available datasets


library(magrittr)

myFiles  %>% nchar %>% order %>% myFiles[.]

# =    4  Base R Pipe  =========================================================

# Since version 4.1, base R now supports a pipe operator without the need
# to load a special package. Such an introductions of external functionality
# into the language is very rare.
#
# Unfortunately it won't (yet) work with the '[' function, so we need to write
# an intermediate function for this example
extract <- function(x, v) {
  return(v[x])
}

myFiles |> nchar() |> order() |> extract(myFiles)


# =    5  Intermediate Assignment  =============================================

# So what's the problem? As you can see, the piped code may be concise and
# expressive. But there is also a large amount of implicit assignment and
# processing going on and that is usually a bad idea because it makes code hard
# to maintain. I am NOT a big fan of the nested syntax, but I don't think that
# replacing it with the pipe makes things much better. My preferred idiom is
# to use intermediate assignments. Only then is it convenient to examine
# the code step by step and validate every single step. And that is the most
# important objective at all: no code is good if it does not compute
# correctly.


x <- nchar(myFiles)
x <- order(x)
myFiles[x]



# =    6  Postscript  ==========================================================

# I tried to write an example that strips all comments from a list of files, and
# another example that finds all files that were not yet updated this year
# (according to the "# Date: in the header). Neither examples can be well
# written without intermediate assignments, or at least sapply() functions
# that are not simpler at all than the intermediate assignment.

# [END]
