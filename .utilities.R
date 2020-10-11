# tocID <- "./.utilities.R"
#
# Miscellaneous R code to suppport the project
#
# Version: 1.4
# Date:    2017-09 - 2020-09
# Author:  Boris Steipe
#
# V 1.4    Maintenance, and new validation utilities
# V 1.3.1  prefix Biostrings:: to subseq()
# V 1.3    load msa support functions
# V 1.2    update database utilities to support 2017 version of JSON sources
# V 1.1    2017 updates for ABC-units
# V 1.0    First code
#
# ToDo:
# Notes:
#
# ==============================================================================


#TOC> ==========================================================================
#TOC> 
#TOC>   Section  Title                                       Line
#TOC> -----------------------------------------------------------
#TOC>   1        SCRIPTS TO SOURCE                             52
#TOC>   2        PACKAGES                                      58
#TOC>   3        DATA & CONSTANTS                              69
#TOC>   4        SUPPORT FUNCTIONS                            116
#TOC>   4.01       objectInfo()                               119
#TOC>   4.02       biCode()                                   147
#TOC>   4.03       sameSpecies()                              181
#TOC>   4.04       validateFA()                               201
#TOC>   4.05       readFASTA()                                309
#TOC>   4.06       writeFASTA()                               344
#TOC>   4.07       pBar()                                     377
#TOC>   4.08       waitTimer()                                399
#TOC>   4.09       fetchMSAmotif()                            427
#TOC>   4.10       H() (Shannon entropy)                      471
#TOC>   4.11       CX() (ChimeraX remote command)             484
#TOC>   5        FUNCTIONS TO CUSTOMIZE ASSIGNMENTS           541
#TOC>   5.01       seal()                                     543
#TOC>   5.02       getMYSPE()                                 547
#TOC>   5.03       selectPDBrep()                             558
#TOC>   5.04       sealKey()                                  593
#TOC>   5.05       selectChi2()                               623
#TOC>   5.06       selectENSP()                               636
#TOC> 
#TOC> ==========================================================================


# =    1  SCRIPTS TO SOURCE  ===================================================

source("./scripts/ABC-dbUtilities.R")
source("./scripts/ABC-writeALN.R")
source("./scripts/ABC-writeMFA.R")

# =    2  PACKAGES  ============================================================

if (! requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest")
}

if (! requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}


# =    3  DATA & CONSTANTS  ====================================================

# cf. https://www.bioinformatics.org/sms/iupac.html
AAVALID  <- "acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY*-"
NUCVALID <- "acgtuACGTU-"
NUCAMBIG <- "acgtACGTryswkmbdhvnRYSWKMBDHVN-"

# A colour palette for amino acid properties
AACOLS <- character()
AACOLS["R"] <- "#5770ff" # Positive
AACOLS["K"] <- "#4785EE" #
AACOLS["H"] <- "#37a1de" #
AACOLS["E"] <- "#ff6f59" # Negative
AACOLS["D"] <- "#ff7391" #
AACOLS["N"] <- "#C9D4FF" # Hydrophilic
AACOLS["Q"] <- "#CADFFC" #
AACOLS["S"] <- "#CBEAF9" #
AACOLS["T"] <- "#CDF5F7" #
AACOLS["Y"] <- "#FBFFC9" # Hydrophobic
AACOLS["W"] <- "#EDFDC8" #
AACOLS["F"] <- "#DFFCC8" #
AACOLS["I"] <- "#D2FBC8" #
AACOLS["L"] <- "#C4FAC7" #
AACOLS["M"] <- "#B7F9C7" #
AACOLS["V"] <- "#A9F8C7" #
AACOLS["A"] <- "#9CF7C7" #
AACOLS["G"] <- "#d2d2d2" # Glycine
AACOLS["C"] <- "#fff963" # Cysteine
AACOLS["P"] <- "#edc06d" # Proline
AACOLS <- gsub("$", "80", AACOLS)  # Make the colors 50% transparent
# barplot(rep(1, 20), col = AACOLS)
# colorRampPalette(c("#fbffc9","#9cf7c7"))(8)

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
                "Wallemia mellicola")


# =    4  SUPPORT FUNCTIONS  ===================================================


# ==   4.01  objectInfo()  =====================================================
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


# ==   4.02  biCode()  =========================================================
biCode <- function(s) {
  # Make a 5 character "biCode" from a binomial name by concatening
  # the uppercased first three letter of the first word and the first
  # two letters of the second word. If there is only one word, we take the
  # first five characters from that. Outputs are padded with "." if necessary.
  # NAs in input are preserved.
  # Parameters:
  #    s   chr  vector of binomial species names
  # Value: chr  vector of biCodes, same length as s, NAs are preserved

  b <- character(length(s))
  s <- gsub("[^a-zA-Z ]", "", as.character(s)) # remove all non-alphabetic
                                               # characters except space
  s <- toupper(s)

  for (i in seq_along(s)) {
    x <- unlist(strsplit(s[i], "\\s+"))
    if (length(x) == 0) { # empty string
      x <- c("", "")
    } else if (length(x) == 1) { # only one string
      x <- c(substr(x, 1, 3), substr(x, 4, 5)) # 3 + 2 with whatever is there
    }
    x <- paste0(x[1:2], "...")  # pad strings

    b[i] <- paste0(substr(x[1], 1, 3), substr(x[2], 1, 2))
  }

  b[is.na(s)] <- NA  # recover NAs from input

  return(b)
}


# ==   4.03  sameSpecies()  ====================================================
sameSpecies <- function(a, b) {
  # Parameters: a, b two vectors that contain
  # binomial species names and maybe additional strain information.
  # Value: a boolean vector, true where the species in a is the same as
  # the species in b.
  # Note: the usual vector recycling applies. Length is not checked.
  a <- gsub("^(\\S+\\s\\S+).*", "\\1", a)
  b <- gsub("^(\\S+\\s\\S+).*", "\\1", b)
  if (any(! grepl("^\\S+\\s\\S+$", a))) {
    stop("\"a\" contains elements that are not binomial names.")
  }
  if (any(! grepl("^\\S+\\s\\S+$", b))) {
    stop("\"b\" contains elements that are not binomial names.")
  }
  return(a == b)
}



# ==   4.04  validateFA()  =====================================================

validateFA <- function(txt) {
  # validates txt according to FASTA assumptions
  # Parameters:
  #   txt   char   a putative vector of FASTA formatted text
  # Value:  invisible(NULL)
  #
  # The function is used for its side-effect of throwing an error of
  # FASTA assumptions are violated in txt.
  # - There may be spacer-lines matching "^\\s*$"
  # - At least one header line
  # - No adjacent header lines
  # - All header lines followed by at least one sequence line
  # - All sequence lines have at least one valid character
  # - Valid characters are AAVALID (which includes valid nucleotides)

  if ( ! any(grepl("^>", txt)) ) {
    stop("no header lines in input")
  }

  isHeader <- grepl("^>", txt)
  isSpacer <- grepl("^\\s*$", txt)
  txtNoSpc <- txt[! isSpacer]

  # adjacent headers
  sel <- isHeader[- length(isHeader)] & isHeader[-1]
  if ( any(sel) ) {
    i <- which(sel)[1]
    stop(sprintf("adjacent header lines in input (lines %d and %d)",
                 i, i+1))
  }

  # invalid character in a line that is not header or spacer
  selA <- isHeader | isSpacer
  selB <- grepl(sprintf("[^%s]", AAVALID), txt)
  if ( any( (! selA) & selB) ) {  # (not header) AND (has invalid character)
    i <- which( (! selA) & selB)[1]
    stop(sprintf("invalid character(s) in sequence (cf. line %d)", i))
  }

  # no spacer should immediately follow a header
  sel <- c(isSpacer, FALSE) & c(FALSE, isHeader)
  if ( any(sel) ) {
    i <- which(sel)[1]
    stop(sprintf("a header has no adjacent sequence (line %d)", i - 1))
  }

  # no header alone on the last line
  if ( which(isHeader)[sum(isHeader)] == length(txt)) {
    stop(sprintf("a header is alone on the last line (line %d)", length(txt)))
  }

  # no spacer in a block of sequence
  # (tested as: not followed by spacer or header)

  if (sum(isSpacer) > 0) {
    sel <- which(isSpacer) + 1
    sel <- sel[sel <= length(txt)]
    if ( ! all(isSpacer[sel] | isHeader[sel]) ) {
      i <- which(! (isSpacer[sel] | isHeader[sel]))[1]
      stop(sprintf("a spacer is followed by sequence (line %d)",
                   which(isSpacer)[i]))
    }

  }

  # all good, if we get to here.
  return(invisible(NULL))
}

if (FALSE) {

  # Should validate
  fa <- c(">abc", "de")
  validateFA(fa)

  fa <- c(">abc", "de", "fgh", ">ojz", "ikl")
  validateFA(fa)

  fa <- c("", ">abc", "de", "fgh", "", ">ojz", "ikl", "")
  validateFA(fa)

  fa <- c(" ", ">abc", "de", "fgh", "   ", ">ojz", "ikl", "\t")
  validateFA(fa)

  # should fail
  fa <- c("abc", "de")
  validateFA(fa)

  fa <- c(">abc", "de", ">fgh", ">ojz", "ikl")
  validateFA(fa)

  fa <- c("", ">abc", "de", "f_h", "", ">ojz", "ikl", "")
  validateFA(fa)

  fa <- c(">abc", " ", "de", "fgh", "   ", ">ojz", "ikl", "\t")
  validateFA(fa)

  fa <- c(" ", ">abc", "de", "fgh", "   ", ">ojz")
  validateFA(fa)

  fa <- c(" ", ">abc", "de", "   ", "fgh", ">ojz", "ikl", "\t")
  validateFA(fa)

}


# ==   4.05  readFASTA()  ======================================================

readFASTA <- function(FA) {
  # Read FASTA formatted text, validate it,
  # return a dataframe of headers and collapsed sequences.
  # Parameters:
  #    FA  chr   Input file name (or text connection)
  # Value:
  #    data.frame
  #       $header char  the FASTA header lines
  #       $ seq   char  the actual sequences
  #
  # Note: if length(FA) is one, it is assumed to be a filename
  #
  # Example:
  #   refAPSES <- readFASTA("./data/refAPSES.mfa")
  #   readFASTA(c("> This", "acdef", "ghi", > That", "k-l"))
  #

  if (length(FA) == 1) { FA <- readLines(FA) }
  validateFA(FA)
  FA <- FA[! grepl("^$", FA)]   # drop all empty lines
  iHead <- grep("^>", FA) # find all headers
  myFA <- data.frame(head = FA[iHead],
                     seq  = character(length(iHead)))

  for (i in seq_along(iHead)) {
    first <- iHead[i] + 1   # first line of each sequence
    last  <- ifelse(i < length(iHead), iHead[i + 1] - 1, length(FA)) # ...last
    myFA$seq[i] <- paste0(FA[first:last], collapse = "")
  }
  return(myFA)
}


# ==   4.06  writeFASTA()  =====================================================
writeFASTA <- function(fa, fn = NULL, width = 60) {
  # Write the contents of dataframe "fa" as a FASTA formatted file.
  # Parameters:
  #    fa      dataframe
  #      $head chr  vector of FASTA headers,
  #      $seq  chr  vector of sequences in one-letter code
  #    fn      chr  filename for output; if NULL (default) the output is
  #                 returned instead
  #    width   int   max number of sequence characters per line of output.
  # Value:
  #           FASTA formatted character vector IF fn was NULL. invisible(NULL)
  #           otherwise.

  out <- character()
  for (i in seq_along(fa$head)) {
    out <- c(out, fa$head[i])                 # add header line
    from <- seq(1, nchar(fa$seq[i]), by = width) # starting indices of chunks
    to <- c((from - 1)[-1], nchar(fa$seq[i]))     # ending indices of chunks
    out <- c(out, substring(fa$seq[i], from, to)) # add chunks to txt
    out <- c(out, "")         # add empty line for better readability
  }
  out <- out[ - length(out)]  # drop the last empty line

  if (length(fn) == 1) {
    writeLines(out, fn)
    return(invisible(NULL))
  } else {
    return(out)
  }
}


# ==   4.07  pBar()  ===========================================================
  pBar <- function(i, l, nCh = 50) {
    # Draw a progress bar in the console
  # i: the current iteration
  # l: the total number of iterations
  # nCh: width of the progress bar
  ticks <- round(seq(1, l-1, length.out = nCh))
  if (i < l) {
    if (any(i == ticks)) {
      p <- which(i == ticks)[1]  # use only first, in case there are ties
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


# ==   4.08  waitTimer()  ======================================================
waitTimer <- function(t, nIntervals = 50) {
  # pause and wait for t seconds and display a progress bar as
  # you are waiting
  t <- as.numeric(t)

  if (t < 0.1) {return(invisible())}

  increment <- t / nIntervals

  bar <- "----:----|"  # One module for the progress bar:
  bar <- rep(bar, ceiling(nIntervals / 10))  # repeat,
  bar <- unlist(strsplit(bar, "")) # split into single characters,
  bar <- bar[1:nIntervals]  # truncate,
  bar <- paste(bar, collapse="") # and collapse.

  cat(sprintf("\nWaiting: |%s\n         |", bar))
  for (i in 1:(nIntervals - 1)) {
    Sys.sleep(increment)
    cat("=")
  }
  Sys.sleep(increment)
  cat("|\n\n")

  return(invisible())
}


# ==   4.09  fetchMSAmotif()  ==================================================
fetchMSAmotif <- function(ali, mot) {
  # Retrieve a subset from ali that spans the sequence in mot.
  # Biostrings package must be installed.
  # Parameters:
  #    ali        MsaAAMultipleAlignment object
  #    mot  chr   substring within ali
  # Value:  AAStringset

  if (class(ali) != "MsaAAMultipleAlignment" &&
      class(ali) != "MsaDNAMultipleAlignment" &&
      class(ali) != "MsaRNAMultipleAlignment") {
    stop("ali has to be an msa multiple alignment object.")
  }

  if (class(mot) != "character") {
    stop("mot has to be a character object.")
  }

  x <- gsub("-", "", as.character(ali))  # pure sequence, no hyphens

  idx <- grep(mot, x)[1] # first sequence containing mot. If no match,
                         # idx becomes NA
  if (is.na(idx)) {
    stop("mot is not a subsequence in ali.")
  }

  # Find the match range
  m <- regexpr(mot, x[idx])
  motifStart <- as.numeric(m)
  motifEnd <- attr(m, "match.length") + motifStart - 1

  # Count characters, skip hyphens ...
  x <- unlist(strsplit(as.character(ali)[idx], ""))
  x <- x != "-"
  x <- as.numeric(x)
  x <- cumsum(x)

  return(Biostrings::subseq(ali@unmasked,
                start = which(x == motifStart)[1], # get the first position
                end   = which(x == motifEnd)[1]))
}


# ==   4.10  H() (Shannon entropy)  ============================================
H <- function(x, N) {
  # calculate the Shannon entropy of the vector x given N possible states
  # (in bits).
  # H(x) = - sum_i(P(x_i) * log2(P(x_i)); 0 * log(0) == 0
  t <- table(x)
  if (missing(N)) { N <- length(t) }
  if (length(t) > N ) { stop("N can't be smaller than observed states.") }
  h <- sum(- (t / length(x)) * log2(t / length(x)))
  return(h)
}


# ==   4.11  CX() (ChimeraX remote command)  ===================================
CX <- function(cmd, port = CXPORT, quietly = FALSE) {
  # send a command to ChimeraX listening on port CXPORT via its REST
  # interface.
  # Parameters:
  #   cmd      char     a ChireaX commandline command
  #   port     int      the portnumber on which ChimeraX is listening
  #   quietly  logical  if FALSE, cat() the contents of the response
  #
  # Value:  the reply by ChimeraX, invisibly.

  CXREST <- sprintf("http://127.0.0.1:%s/run?", CXPORT)

  cmd <- gsub("(^\\s+)|(\\s+$)", "", cmd)  # trim whitespace
  # percent encode reserved characters
  cmd <- gsub("%",   "%25", cmd)          #   %
  cmd <- gsub("#",   "%23", cmd)          #   #
  cmd <- gsub("/",   "%2F", cmd)          #   /
  cmd <- gsub(":",   "%3A", cmd)          #   :
  cmd <- gsub("@",   "%40", cmd)          #   @
  cmd <- gsub(",",   "%2C", cmd)          #   ,
  cmd <- gsub("\\*", "%2A", cmd)          #   *
  cmd <- gsub("\\?", "%3F", cmd)          #   ?
  cmd <- gsub("!",   "%21", cmd)          #   !
  cmd <- gsub("=",   "%3D", cmd)          #   =
  cmd <- gsub("\\(", "%28", cmd)          #   (
  cmd <- gsub("\\)", "%29", cmd)          #   )
  cmd <- gsub("\\[", "%5B", cmd)          #   [
  cmd <- gsub("\\]", "%5D", cmd)          #   ]
  cmd <- gsub("&",   "%26", cmd)          #   &
  cmd <- gsub("\\+", "%2B", cmd)          #   +

  cmd <- gsub("\\s+", "+", cmd)            # whitespace to "+"
  cmd <- URLencode(cmd)                    # encode special characters
  cmd <- paste0(CXREST, "command=", cmd, collapse = "")

  r <- httr::GET(cmd)

  if (! r$status_code == 200) {
    stop("ChimeraX returned status code %d", r$status_code)
  }

  if (length(r$content) == 0) {
    reply <- ""
  } else {
    reply <- rawToChar(r$content)
  }

  if (quietly == FALSE) {
    cat(reply)
  }

  return(invisible(reply))

}


# =    5  FUNCTIONS TO CUSTOMIZE ASSIGNMENTS  ==================================

# ==   5.01  seal()  ===========================================================
seal <- function(x.1L) { .Call(digest:::digest_impl,x.1L,3L,-1L,-0,-0,-0) }


# ==   5.02  getMYSPE()  =======================================================
getMYSPE <- function(x) {
  dat <- readRDS("./data/sDat.rds")
  map <- readRDS("./data/MYSPEmap.rds")
  key <- gsub(".+(....).$","\\1", x)
  x <- dat$species[map[key,"iMYSPE"]]
  names(x) <- dat$taxID[map[key,"iMYSPE"]]
  return(x)
}


# ==   5.03  selectPDBrep()  ===================================================
selectPDBrep <- function(n, forCredit = FALSE) {
  # Select n PDB IDs from a list of high-resolution, non-homologous, single
  # domain, single chain structure files that represent a CATH topology
  # group.
  # Parameters:
  #   n     num     number of IDs to return
  #   seed  num     a seed for the RNG
  #
  # Value:          char  PDB IDs
  #
  # Note: the list is loaded from an .rds file in the "./data" directory.

  if (forCredit) {
    seed <- myStudentNumber
  } else {
    seed <- as.integer(Sys.time())
    cat("NOTE: This selection will not validate for a course submission.\n")
    cat("      If you intend to use it for an assignment task, invoke\n")
    cat("      this function like \"selectPDBrep(n, forCredit = TRUE)\".\n\n")
  }

  pdbRep <- readRDS("./data/pdbRep.rds")  # loads pdbRep

  if (n > length(pdbRep)) {
    stop(sprintf("There are only %d PDB IDs in the table to choose from.",
                 length(pdbRep)))
  }
  oldSeed <- .Random.seed
  set.seed(seed)
  PDBset <- sample(pdbRep, n)
  .Random.seed <- oldSeed
  return(PDBset)
}

# ==   5.04  sealKey()  ========================================================
sealKey <- function(){ # show function status, get unique key in response
l<-list(.Random.seed,Sys.time,sys.frame,sys.nframe,sys.parent,capture.output,
str,as.integer,gsub,get,11,function(s,q=1,k=l[[11]]){
s<-as.integer(charToRaw(s))-32;set.seed(k);y<-sample(0:94,length(s),r=T)*q;
.Random.seed<-l[[1]];s[s<F|s>(95-T)]<-F;s<-(s+y) %%(9^2+3^2+2^2);s<-s+32;
s<-rawToChar(as.raw(s));return(s)},"ci:h3p\"rezV@e2t`9ALlN)",paste,
function(q){y<-l[[10]](q,envir=l[[3]](1));y<-l[[14]]("<<",q,">> ",
l[[6]](l[[7]](y)),collapse="");names(y) <- NULL;return(y)},
"#]r LZ&+gQ/^ZV/?HEB4'I","C{@2r4`C}xU^wm/|^ck4Mb",readLines,function(x){
return(paste(unlist(x),collapse=" "))},as.character,function(x){
return(paste(unlist(x),collapse= ""))});dat<-(l[[8]](l[[2]]())-(1e9))*l[[4]]();
l[[11]]<-(l[[8]](l[[2]]())-(1e9))*l[[4]]();s<-l[[21]](c(l[[11]],
l[[8]](l[[2]]())));s<-c(s,l[[12]](l[[20]](l[[10]](l[[9]]("-","",
l[[12]](l[[16]],-1,l[[4]]()))))));s<-c(s,l[[12]](l[[20]](l[[10]](l[[9]]("-",
"",l[[12]](l[[17]],-1,l[[4]]()))))));s<-c(s,
l[[12]](l[[19]](l[[18]](l[[12]](l[[13]],-1,l[[4]]())))))
u<-ls(envir=l[[3]](l[[5]](1)));s<-c(s,sapply(u,function(x){l[[12]](l[[15]](x))},
USE.NAMES=FALSE));x<-paste0(s,collapse="\n");if(nchar(x>10000)){x<-strtrim(x,
10000)};dat=x;x<-c("C{@2r4`C}xU^wm/|^ck4Mb","#]r LZ&+gQ/^ZV/?HEB4'I")
# ==============================================================================
response <- httr::POST("http://steipe.biochemistry.utoronto.ca/abc/seal.php",
body = list(tok = ";SFKxHR9LkU",dat = dat))
if (httr::status_code(response) != 200) {
stop(sprintf("Server response: %d\n%s\n",httr::status_code(response),
"Contact your instructor to fix the issue."))}
lines <- unlist(strsplit(httr::content(response, "text"), "\\n"))
print(lines[grep("sealKey", lines)]);return(invisible(NULL))}


# ==   5.05  selectChi2()  =====================================================
selectChi2 <- function() {
  # Select one random Amino acid from those that have a Chi2 angle

  oldSeed <- .Random.seed
  set.seed(myStudentNumber)
  AA <- sample(c("Asp", "Glu", "Phe", "His", "Ile", "Lys", "Leu",
                 "Met", "Asn", "Gln","Arg", "Trp", "Tyr"))
  .Random.seed <- oldSeed
  cat(sprintf("  Chi1/Ch2: Use \"%s\".  <%s>\n", AA[4], seal(AA)))
}


# ==   5.06  selectENSP()  =====================================================
selectENSP <- function(x) {
  oldSeed <- .Random.seed
  set.seed(myStudentNumber)
  x <- sample(x[order(x)])
  .Random.seed <- oldSeed
  cat(sprintf("seal: %s\n", seal(paste0(x,collapse=""))))
  return(x)
}



# [END]
