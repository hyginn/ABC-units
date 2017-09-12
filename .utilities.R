# .utilities.R
#
# Miscellaneous R code to suppport the project
#
# Version: 1.0
# Date:    2016 12
# Author:  Boris Steipe
#
# V 1.0    First code
#
# ToDo:
# Notes:
#
# ==============================================================================

objectInfo <- function(x) {
    # Function to combine various information items about R objects
    #
    # Input: an R object
    # Value: none - prints information as side-effect

    cat("object contents:")
    print(x, digits = 22)  # print value at maximal precision

    cat("\nstructure of object:\n")
    str(x)

    if (! is.list(x)) { # Don't use cat() if x is a list. cat() can't handle lists.
        cat("\nmode:   ", mode(x), "\n")
        cat("typeof: ", typeof(x), "\n")
        cat("class:  ", class(x), "\n")
    }

    # if the object has attributes, print them too
    if (! is.null(attributes(x))) {
        cat("\nattributes:\n")
        attributes(x)
    }
    # Done
}

# ====== Utilities

biCode <- function(s) {
  # make a 5 character code from a binomial name by concatening
  # the first three letter of the first word and the first
  # two letters of the second word.
  #
  # There are many ways to do this, here we assemble the two parts
  # in a loop, this way the function is vectorized and can
  # work on a large vector of names.
  b <- character()
  for (i in 1:length(s)) {
    b[i] <- sprintf("%s%s",
                    toupper(unlist(substr(s[i], 1, 3))),
                    toupper(unlist(substr(strsplit(s[i], "\\s+")[[1]][2],
                                          1, 2))))
  }
  return(b)
}


pBar <- function(i, l, nCh = 50) {
  # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)
      p1 <- paste(rep("#", p), collapse = "")
      p2 <- paste(rep("-", nCh - p), collapse = "")
      cat(sprintf("\r|%s%s|", p1, p2))
      flush.console()
    }
  }
  else { # done
    cat("\n")
  }
}

# ====== DATA ==================================================================


# 10 species of fungi for reference analysis.
# http://steipe.biochemistry.utoronto.ca/abc/index.php/Reference_species_for_fungi
REFspecies <- c("Aspergillus nidulans",
                "Bipolaris oryzae",
                "Coprinopsis cinerea",
                "Cryptococcus neoformans",
                "Neurospora crassa",
                "Puccinia graminis",
                "Saccharomyces cerevisiae",
                "Schizosaccharomyces pombe",
                "Ustilago maydis",
                "Wallemia mellicola"
)


# [END]
