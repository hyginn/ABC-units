# .init.R
# Functions to initialize this collection of learning units
# Boris Steipe
# ====================================================================

# Create a local copy of myScript.R if required, and not been done yet.
if (! file.exists("myScript.R") && file.exists(".tmp.R")) {
    file.copy(".tmp.R", "myScript.R")
}

source(".utilities.R")

file.edit("ABC-units.R")

# [End]
