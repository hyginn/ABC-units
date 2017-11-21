# ABC-writeALN.R
#
# ToDo:    calculate consensus line
#          append sequence numbers
# Notes:
#
# ==============================================================================


writeALN <- function(ali,
                     range,
                     note = "",
                     myCon = stdout(),
                     blockWidth = 60) {
  # Purpose:
  #     Write an MsaAAMultipleAlignment or AAStringSet object to stdout() or
  #     a file in multi-FASTA format.
  # Version: 2.0
  # Date:    2017 10
  # Author:  Boris Steipe
  #
  # Parameters:
  #     ali             MsaAAMultipleAlignment or AAStringSet or character
  #                       vector.
  #     range      num  a two-integer vector of start and end positions if
  #                       only a range of the MSA should be written, e.g.
  #                       a domain. Defaults to the full alignment length.
  #     note       chr  a vector of character that is appended to the name
  #                       of a sequence in the FASTA header. Recycling of
  #                       shorter vectors applies, thus a vector of length one
  #                       is added to all headers.
  #     myCon           a connection (cf. the con argument for writeLines).
  #                       Defaults to stdout()
  #     blockWidth int  width of sequence block. Default 80 characters.
  # Value:
  #     NA   the function is invoked for its side effect of printing an
  #          alignment to stdout() or file.

  blockWidth <- as.integer(blockWidth)
  if (is.na(blockWidth)) {
    stop("PANIC: parameter \"blockWidth\" must be numeric.")
  }
  if (blockWidth < 1) {
    stop("PANIC: parameter \"blockWidth\" must be greater than zero.")
  }
  if (blockWidth > 60) {
    warning("Programs that read CLUSTAL format might not expect blockWidth > 60.")
  }

  # Extract the raw data from the objects depending on their respective class
  # and put it into a named vector of strings.

  # Extract XStringSet from MsaXMultipleAlignment ...
  if (class(ali) == "MsaAAMultipleAlignment" |
      class(ali) == "MsaDNAMultipleAlignment" |
      class(ali) == "MsaRNAMultipleAlignment") {
      ali <- ali@unmasked
  }

  # Process XStringSet
  if (class(ali) == "AAStringSet" |
      class(ali) == "DNAStringSet" |
      class(ali) == "RNAStringSet") {
    sSet <- as.character(ali) # we use as.character(), not toString() thus
                              # we don't _have_ to load Biostrings
  } else if (class(ali) == "character") {
    sSet <- ali
  } else {
    stop(paste("Input object of class",
               class(ali),
               "can't be handled by this function."))
  }

  if (missing(range)) {
    range <- 1
    range[2] <- max(nchar(sSet))
  } else {
    range <- as.integer(range)
    if(length(range) != 2 ||
       any(is.na(range)) ||
       range[1] > range[2] ||
       range[1] < 1) {
      stop("PANIC: \"range\" parameter must contain valid start and end index.")
    }
  }

  # Right-pad any sequence with "-" that is shorter than ranges[2]
    for (i in seq_along(sSet)) {
      if (nchar(sSet[i]) < range[2]) {
        sSet[i] <- paste0(sSet[i],
                          paste0(rep("-", range[2] - nchar(sSet[i])),
                                 collapse = ""))
      }
    }

  # Right-pad sequence names
  sNames <- names(sSet)
  len <- max(nchar(sNames)) + 2 # longest name plus two spaces
  for (i in seq_along(sNames)) {
    sNames[i] <- paste0(sNames[i],
                      paste0(rep(" ", len - nchar(sNames[i])),
                             collapse = ""))
  }


  # Process each sequence
  txt <- paste0("CLUSTAL W format. ", note)
  txt[2] <- ""

  iStarts <- seq(range[1], range[2], by = blockWidth)
  iEnds <- c((iStarts[-1] - 1), range[2])

  for (i in seq_along(iStarts)) {
    for (j in seq_along(sSet)) {
      txt <- c(txt,
               paste0(sNames[j], substring(sSet[j], iStarts[i], iEnds[i])))
    }
    txt <- c(txt, "")  # append a blank consenus line
    txt <- c(txt, "")  # append a separator line
  }

  writeLines(txt, con= myCon)

}

# ====  TESTS  =================================================================
# Enter your function tests here...

if (FALSE) {
  # test ...
}



# [END]
