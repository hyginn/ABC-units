# tocID <- "ABC-Install_all_packages.R"
#
# Purpose:  A Bioinformatics Course:
#              Installing all packages in this course
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
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                          Line
#TOC> ----------------------------------------------
#TOC>   1        Packages                         33
#TOC>   2        CRAN packages                    98
#TOC>   3        Bioconductor packages           127
#TOC>   4        Other package sources           142
#TOC>   5        Updating packages               148
#TOC>
#TOC> ==========================================================================


# =    1  Packages  ============================================================

# Much of R's functionality is contributed in packages: bundles of R scripts
# or code in other languages, pre-configured objects, and datasets. Making this
# functionality available is often done by issuing a library(<package-name>)
# command, however this is not the preferred way, since it may override other
# R functions and it makes it harder to understand where the source code of
# a particular function is located. In this course we call the function name
# prefixed with the package name and two colons:
#   <package-name>::<function-name>()
# This is the preferred way, since it is explicit.
#
# Regardless of which idiom one uses to call the actual function, the package
#  needs to be "installed" first, i.e. the code must have been downloaded
# from CRAN, or using the BiocManager::install() function.
#
# This script contains download commands for all packages that are used in the
# course. You can execute the script line by line (or even source the entire
# script) to make sure all packages can be installed on your computer. Just
# one reminder: if you are ever asked to install from source, the correct
# answer is usually "no" - except if you really know what you are doing and why.
#
# Once packages are installed you can get additional information about
# the contents of a package with the commands:
#  library(help=<package-name>)       # basic information
#  browseVignettes("<package-name>")  # available vignettes
#  data(package = "<package-name>")   # available datasets
#
#  ... and you can load data sets with:
#  data(<data-set-name>, package = "<package-name>")
#
#  All packages here are installed only when they have not been installed
#  before, using the following idiom:
#
#     if (! requireNamespace("<package-name>", quietly=TRUE)) {
#       install.packages("<package-name>")
#     }
#
#  ... or its BiocManager::install() equivalent:
#
# if (! requireNamespace("<bioconductor-package-name>", quietly=TRUE)) {
#   BiocManager::install("<bioconductor-package-name>")
# }
#
#  If you want to _force_ a re-installation of the package, simply issue
#  the install.packages("<package-name>") command on its own. For compactness
#  we wrap the idiom into a function, which can also switch between CRAN
#  and BIOconductor sources:

installIfNeeded <- function(package, s = "CRAN") {
  # s: "CRAN" or "BIO"
  if (s == "CRAN") {
    if (! requireNamespace(package, quietly=TRUE)) {
      install.packages(package)
    }
  } else if (s == "BIO") {
    if (! requireNamespace("BiocManager", quietly=TRUE)) {
      install.packages("BiocManager")
    }
    if (! requireNamespace(package, quietly=TRUE)) {
      BiocManager::install(package)
    }
  } else {
    stop(sprintf("Unknown source \"%s\".", s))
  }
}


# =    2  CRAN packages  =======================================================

installIfNeeded("ape")
installIfNeeded("BiocManager")
installIfNeeded("bio3d")
installIfNeeded("evd")
installIfNeeded("ggseqlogo")
installIfNeeded("ggtern")
installIfNeeded("hexbin")
installIfNeeded("httr")
installIfNeeded("igraph")
installIfNeeded("jsonlite")
installIfNeeded("magrittr")
installIfNeeded("MASS")
installIfNeeded("microbenchmark")
installIfNeeded("phangorn")
installIfNeeded("plotly")
installIfNeeded("plotrix")
installIfNeeded("profvis")
installIfNeeded("robustbase")
installIfNeeded("RColorBrewer")
installIfNeeded("Rphylip")
installIfNeeded("rvest")
installIfNeeded("seqinr")
installIfNeeded("stringi")
installIfNeeded("taxize")
installIfNeeded("testthat")
installIfNeeded("xml2")

# =    3  Bioconductor packages  ===============================================

installIfNeeded("Biobase",       s = "BIO")
installIfNeeded("biomaRt",       s = "BIO")
installIfNeeded("Biostrings",    s = "BIO")
installIfNeeded("DECIPHER",      s = "BIO")
installIfNeeded("GEOquery",      s = "BIO")
installIfNeeded("GOSim",         s = "BIO")
installIfNeeded("limma",         s = "BIO")
installIfNeeded("msa",           s = "BIO")
installIfNeeded("org.Sc.sgd.db", s = "BIO")
installIfNeeded("prada",         s = "BIO")
installIfNeeded("topGO",         s = "BIO")


# =    4  Other package sources  ===============================================

# Using sources other than CRAN or Bioconductor to download general-purpose
# programs that run on your computer is not generally recommended.


# =    5  Updating packages  ===================================================

# From time to time, update CRAN packages with the following command ...

update.packages()

# ... and also update Bioconductor packages as follows:

BiocManager::install()

# [END]
