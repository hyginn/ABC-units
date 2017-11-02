# ABC-writeMFA.R
#
# ToDo:
# Notes:  2.1  bugfix: empty notes caused superfluous blank after header.
#
#
# ==============================================================================


writeMFA <- function(ali,
                     range,
                     note = "",
                     myCon = stdout(),
                     blockWidth = 80) {
  # Purpose:
  #     Write an MsaAAMultipleAlignment or AAStringSet object to stdout() or
  #     a file in multi-FASTA format.
  # Version: 2.1
  # Date:    2017  10
  # Author:  Boris Steipe
  #
  # Parameters:
  #     ali             MsaAAMultipleAlignment or AAStringSet or character
  #                       vector
  #     range      num  a two-integer vector of start and end positions if
  #                       only a range of the MSA should be written, e.g.
  #                       a domain. Defaults to the full sequence length.
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
  if (blockWidth < 1){
    stop("PANIC: parameter \"blockWidth\" must be greater than zero.")
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

  # Process each sequence
  txt <- character()
  if (note != "") {  # construct header line
    headers <- paste(names(sSet), note)
  } else {
    headers <- names(sSet)
  }

  for (i in seq_along(sSet)) {

    # output FASTA header
    txt <- c(txt, sprintf(">%s", headers[i]))

    # output the sequence in blocks of blockWidth per line ...
    iStarts <- seq(range[1], range[2], by = blockWidth)
    iEnds <- c((iStarts[-1] - 1), range[2])

    thisSeq <- substring(sSet[i], iStarts, iEnds)  # collect all blocks
    thisSeq <- thisSeq[! nchar(thisSeq) == 0]      # drop empty blocks
    txt <- c(txt, thisSeq)

    txt <- c(txt, "")  # append an empty line for readability
  }

  writeLines(txt, con= myCon)

}

# ====  TESTS  =================================================================
# Enter your function tests here...

if (FALSE) {
  # test ...
}



# [END]
