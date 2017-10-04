# .init.R
# Functions to initialize this collection of learning units
# Boris Steipe
# ====================================================================

# Create a local copy of myScript.R if required, and not been done yet.
if (! file.exists("myScript.R") && file.exists(".tmp.R")) {
    file.copy(".tmp.R", "myScript.R")
}

# If it doesn't exist yet, set up a profile:
if (! file.exists(".myProfile.R")) {
  # setup profile data
  cat("\nPlease enter the requested values correctly, no spaces, and\n")
  cat("press <enter>.\n")
  e <- readline("Please enter your UofT eMail address: ")
  n <- readline("Please enter your Student Number: ")

  conn <- file(".myProfile.R")
  writeLines(c(sprintf("myEMail <- \"%s\"", e),
               sprintf("myStudentNumber <- %d", as.numeric(n))),
             conn)
  close(conn)
  rm(e, n, conn)
}

# Patch YFO -> MYSPE if necessary:
tmp <- readLines(".myProfile.R")
if (length(grep("^YFO", tmp)) > 0) {
  idx <- grep("^YFO", tmp)
  tmp[idx] <- gsub("^YFO", "MYSPE", tmp[idx])
  writeLines(tmp, ".myProfile.R")
}
rm(tmp)

source(".myProfile.R")

source(".utilities.R")

file.edit("ABC-units.R")

# [End]
