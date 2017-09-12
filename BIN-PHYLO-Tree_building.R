# ___ID___ .R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the ___ID___ unit.
#
# Version:  0.1
#
# Date:     2017  08  28
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.1    First code copied from 2016 material.

#
# TODO:
#
#
# == DO NOT SIMPLY  source()  THIS FILE! =======================================

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
# going on. That's not how it works ...

# ==============================================================================

# = 1 ___Section___

# ==============================================================================
#        PART FOUR: Calculating trees
# ==============================================================================

# Follow the instructions found at phylip's home on the Web to install. If you
# are on a Windows computer, take note of the installation directory.

# After you have installed Phylip on your computer, install the R package that
# provides an interface to the Phylip functions.

if (!require(Rphylip, quietly=TRUE)) {
  install.packages("Rphylip")
  library(Rphylip)
}

# This will install RPhylip, as well as its dependency, the package "ape".

# The next part may be tricky. You will need to figure out where
# on your computer Phylip has been installed and define the path
# to the proml program that calculates a maximum-likelihood tree.

# On the Mac, the standard installation places a phylip folder
# in the /Applications directory. That folder contains all the
# individual phylip programs as <name>.app files. These are not
# the actual executables, but "app" files are actually directories
# that contain the required resources for a program to run.

# The executable is in a subdirectory and you can point Rphylip
# directly to that subdirectory to find the program it needs:
# PROMLPATH <- "/Applications/phylip-3.695/exe/proml.app/Contents/MacOS"

# On Windows you need to know where the rograms have been installed, and you
# need to specify a path that is correct for the Windows OS. Find the folder
# that is named "exe", and right-click to inspect its properties. The path
# should be listed among them.

# If the path looks like "C:\Users\Meng\Programs\phylip-3.695\exe", then your
# assignment has to be
# PROMLPATH <- "C:/Users/Meng/Programs/phylip-3.695/exe"
# (Note: "/", not "\")

# I have heard that your path must not contain spaces, and it is prudent to
# avoid other special characters as well.

# If you are running Linux I trust you know what to do. It's probably
# something like
# PROMLPATH <- "/usr/local/phylip-3.695/bin"

# Confirm that the settings are right.
PROMLPATH                # returns the path
list.dirs(PROMLPATH)     # returns the directories in that path
list.files(PROMLPATH)    # lists the files

# If "proml" is NOT among the files that the last command returns, you
# can't continue.

# Now read the mfa file you have saved, as a "proseq" object with the
# read.protein() function of the RPhylip package:

apsIn <- read.protein("APSES.mfa")
apsIn <- read.protein("~/Desktop/APSES_HISCA.mfa")

# ... and you are ready to build a tree.

# Building maximum-likelihood trees can eat as much computer time
# as you can throw at it. Calculating a tree of 48 APSES domains
# with default parameters of Rproml() runs for more than half a day
# on my computer. But we have only twelve sequences here, so the
# process will take us about 5 to 10 minutes.

apsTree <- Rproml(apsIn, path=PROMLPATH)



# A quick first look:

plot(apsTree)



# = 1 Tasks




# [END]
