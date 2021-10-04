# .Rprofile
#
# This script is automatically executed on startup
# ==============================================================================

init <- function() {

  # Create a local copy of myScript.R if not done yet.
  if (! file.exists("myScript.R") && file.exists(".tmp.R")) {
    file.copy(".tmp.R", "myScript.R")
    cat("A new file \"myScript.R\" was created. You can use it for\n")
    cat("notes and code experiments.\n\n")
  }

  cat("\n\n")
  cat("Please open the file \".myProfile.R\" (click on the file-name in the\n")
  cat("\"files\" pane), edit it and save it.\n")
  cat("Then click the checkbox, and use the More -> Move... dialogue\n")
  cat("to move it into the \"myScripts\" folder.\n\n")

  file.edit("ABC-units.R")
  return(invisible(NULL))
}

if (! file.exists("./myScripts/.myProfile.R")) {
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
  cat("    =================")
  cat("\n\n")
  cat("        WELCOME !\n")
  cat("\n")
  cat("  Type  'init()'  to begin\n\n")
  cat("\n")
  cat("    =================")
  cat("\n\n")

} else {  # local profile exists ... validate state:
  cat("\n\nLoading local functions ...")

  source(".utilities.R")  # local profile appears sane, source utilities
  source("./myScripts/.myProfile.R")

  if (! exists("myEMail")) {  # ... has eMail been defined?
    cat("ERROR !\n")
    cat("=======\n")
    cat("The file \"./myScripts/.myProfile.R\" exists, but\n")
    cat("the variable \"myEMail\" was not loaded.\n")
    cat("Please contact your instructor to continue.\n\n")
  }
  if (! exists("myStudentNumber")) {  # ... has the Student Number been defined?
    cat("ERROR !\n")
    cat("=======\n")
    cat("The file \"./myScripts/.myProfile.R\" exists, but\n")
    cat("the variable \"myStudentNumber\" was not loaded.\n")
    cat("Please contact your instructor to continue.\n\n")
  }
  if (! grepl("^(100.{7})|(99.{7})$", as.character(myStudentNumber))) {
    cat("ERROR !\n")                 # is the Student Number valid?
    cat("=======\n")
    cat("The file \"./myScripts/.myProfile.R\" exists, but\n")
    cat("your Student Number could not be validated.\n")
    cat("Please examine the file \"./myScripts/.myProfile.R\"\n")
    cat(" and fix the problem or contact your instructor to continue.\n\n")
  }


  if (! exists("MYSPE")) {  # if MYSPE has not yet been defined, define it now
                            # ... and write it into the profile.
       prf <- readLines("./myScripts/.myProfile.R")
       iEmail <- grep("^\\s*myStudentNumber\\s*<-", prf)
       out <- prf[1:iEmail]
       out <- c(out, sprintf("MYSPE <- \"%s\" ",
                             getMYSPE(myStudentNumber)))
       out <- c(out, prf[(iEmail+1):length(prf)])
       writeLines(out, "./myScripts/.myProfile.R")

       cat("\n")
       cat(sprintf("MYSPE (%s) was added to \"./myScripts/.myProfile.R\"\n\n",
                   getMYSPE(myStudentNumber)))
       MYSPE <- getMYSPE(myStudentNumber)  # ... define it for immediate use
       rm(prf, iEmail, out)                # cleanup
  }
  cat("... done.\n\n")
}

if (default.stringsAsFactors()) {
  cat("WARNING.\n")
  cat("========\n")
  cat("Your default \"stringsAsFactors\" parameter is set to \"TRUE\".\n")
  cat("This will break some of the code.\n")
  cat("Please contact your instructor to troubleshoot and fix this issue.\n")
  cat("\n")
}

errText <- list()
errText[["noProfileFile"]] <- '
Your PROFILE FILE does not exist. This problem must be fixed to continue.

  The code expects the file "./myScripts/.myProfile.R" to exist and to
  contain your correct eMail address and student number. Detailed
  instructions were given when you first ran the init() command.

  Try running init() again and follow the instructions. Reload youR RStudio
  session and start over with this file.

  If this does not fix the problem, ask for help.
'

errText[["noStudentNumber"]] <- '
Your STUDENT NUMBER has not been defined. This problem must be fixed to continue.

  The code expects the file "./myScripts/.myProfile.R" to exist and to
  contain your correct eMail address and student number. This file gets
  sourced when you start a new R-session, but since you see this error
  message there was a problem.

  Perhaps you need to restart your R-session. Try closing the RStudio
  project and reopening it from the File > Recent Projects menu.

  Perhaps there was a syntax error in your file. Then not all the
  instructions in the file are executed. Check the file: is your
  email perhpas not defined? Or did you type it without qwuoataion
  marks?

  Try fixing problems, and then restart R as described above.

  If none of this fixes the problem, ask for help.
'

# [END]
