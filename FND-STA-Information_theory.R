# tocID <- "FND-STA-Information_theory.R"
#
# ==============================================================================
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the FND-STA-Information_theory unit.
#
# Version:  0.2.1
#
# Date:     2017 - 2021
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           0.2.1  Maintenance
#           0.2    Under development
#           0.1    First code copied from 2016 material.
#
#
# TODO:
#
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
#TOC>   Section  Title                  Line
#TOC> --------------------------------------
#TOC>   1        ___Section___            39
#TOC> 
#TOC> ==========================================================================


# =    1  ___Section___  =======================================================

# What level of information is "significant"

# Assume the background distribution is the database frequencies of
# amino acids:

AAref <- numeric()  # Uniprot frequencies October 2017, slightly adjusted to
# sum to 1.0
AAref["A"] <- 0.0904
AAref["C"] <- 0.0123
AAref["D"] <- 0.0545
AAref["E"] <- 0.0617
AAref["F"] <- 0.0394
AAref["G"] <- 0.0724
AAref["H"] <- 0.0221
AAref["I"] <- 0.0573
AAref["K"] <- 0.0504
AAref["L"] <- 0.0986
AAref["M"] <- 0.0240
AAref["N"] <- 0.0392
AAref["P"] <- 0.0486
AAref["Q"] <- 0.0381
AAref["R"] <- 0.0570
AAref["S"] <- 0.0673
AAref["T"] <- 0.0558
AAref["V"] <- 0.0686
AAref["W"] <- 0.0129
AAref["Y"] <- 0.0294
sum(AAref)

# Function to calculate Shannon entropy
H <- function(pmf) {
  # Calculate Shannon entropy
  # Parameters:
  #   pmf (numeric) probability mass function: a vector of states and
  #                 associated probabilities. Each element of
  #                 pmf must be in (0, 1] and sum(pmf) must be 1.
  # Value:
  #   Shannon entropy in bits.
  # Examples:
  #   H(c(A=0.25, C=0.25, G=0.25, T=0.25))  # 2 bits entropy in a random
  #                                         # nucleotide sequence
  #   H(1)     # If all elements are the same, entropy is zero
  #
  if (any(pmf <= 0 | pmf > 1) || isFALSE(all.equal(1.0, sum(pmf)))) {
    stop("Input is not a discrete probability distribution.")
  }
  H <- -sum(pmf * (log(pmf) / log(2)))
  return(H)
}

# Why use all.equal()? Exact comparisons with floating point numbers are
# brittle. Consider for example:
1/6 + 1/6 + 1/6 + 1/6 + 1/6 + 1/6 == 1
print(1/6 + 1/6 + 1/6 + 1/6 + 1/6 + 1/6, digits = 22) # 0.9999999999999998889777
# all.equal() tests for _near_ equality with tolerance of ~ 1.5e-8



# Entropy of the database frequencies (in bits):
(Href <- H(AAref))

# for comparison: entropy if all amino acids are equiprobable
H(rep(0.05, 20))


# Set up a simulation to estimate the distribution of Information values
# from random sequences drawn from AAref. This is the distribution for the
# statistical null hypothesis:
nObs <- 15                      # number of observations (e.g aligned sequences)
# nObs <- 80
nTrials <- 10000                # number of trials
IObs <- numeric(nTrials)        # vector to store Information in each trial
simCounts <- numeric(20)        # vector to tabulate our information ...
names(simCounts) <- names(AAref)# ... with the names of AAref


for (i in 1:nTrials) {  # simulate ...

  # sample AAref letters, nObs times, with the probabilities of AAref:
  AAobs <- sample(names(AAref), size = nObs, prob = AAref, replace = TRUE)

  x <- table(AAobs)                            # table simulated observations
  simCounts[1:20] <- rep(0, length(simCounts)) # initialize simCounts to 0
  simCounts[names(x)] <- x                     # overwrite with observed counts
  simCounts <- simCounts + 0.5                 # add Jeffreys' pseudocounts
  Hobs <- H(simCounts/sum(simCounts))          # counts to frequency, calc. H
  IObs[i] <- Href - Hobs                       # store information
}

# evaluate
hist(IObs, col = "#C9F4E3", xlim = c(-0.2, 1.0), breaks = 25)
abline(v = quantile(IObs, c(0.05, 0.95)), col = "#AA00CC")

# The purple lines are drawn at the 5% quantiles of the Iobs distributions -
# i.e. an actual observation that lies outside the purple lines is deemed
# "significant"(1)(2). Of course, this is only true to the degree that the
# database frequencies are a valid model for the null-hypothesis on the
# sequence position we are considering here.

#  (1) If we use 5% quantiles, this means a value is significantly larger
#      than expected, and we ignore cases when the value is < 0; if we
#      consider both smaller and larger values, we need to use 2.5% quantiles,
#      since 5% of all observations lie outside the 0.025 and 0.975
#      quantiles.
#
#  (2) For an actual observation of counts, we calculate its observed
#      _empirical_p_Value_ as (nCounts + 1)/(nTotal + 1).


# You can probably now appreciate that information is a bit of a shortcut for
# biological sequences, and does not really take the different inherent
# frequencies based on the character of the amino acids into account. For
# example, L is the most frequent and C is the least frequent, but if we have an
# alignment of 1000 sequences and we see that the frequencies for L and C are
# swapped, that would be _very_ surprising - nevertheless, the information would
# be 0. In order to take that into account, we should actually compute
# Kullback-Leibler divergences.


# Swap C and L frequencies
p <- AAref
q <- AAref
q["L"] <- AAref["C"]
q["C"] <- AAref["L"]
H(p)
H(q)

KLdiv <- function(p, q) {
  # p and q are two pmfs of discrete probability distributions
  # with the same outcomes, which are nowhere 0.
  # Value:  Kullback-Leibler divergence  sum(p * log( p / q))).

  if (length(p) != length(q)) {
    stop("PANIC: input vector lengths differ!")
  }
  if (any(c((p == 0), (q == 0)))) {
    stop("PANIC: 0's found in input vectors!")
  }

  return(sum(p * log( p / q )))
}

KLdiv(p, p)
KLdiv(p, q)


nObs <- 15                      # number of observations (e.g aligned sequences)
# nObs <- 80
nTrials <- 10000                # number of trials
KLdivObs <- numeric(nTrials)        # vector to store Information in each trial
simCounts <- numeric(20)        # vector to tabulate our information ...
names(simCounts) <- names(AAref)# ... with the names of AAref


for (i in 1:nTrials) {  # simulate ...

  # sample AAref letters, nObs times, with the probabilities of AAref:
  AAobs <- sample(names(AAref), size = nObs, prob = AAref, replace = TRUE)

  x <- table(AAobs)                            # table simulated observations
  simCounts[1:20] <- rep(0, length(simCounts)) # initialize simCounts to 0
  simCounts[names(x)] <- x                     # overwrite with observed counts
  simCounts <- simCounts + 0.5                 # add Jeffreys' pseudocounts
  simCounts <- simCounts/sum(simCounts)        # counts to frequency
  KLdivObs[i] <- sum(simCounts * log( simCounts / AAref )) # store KLdiv
}

# evaluate
hist(KLdivObs, col = "#C9F4E3", breaks = 25)
abline(v = quantile(KLdivObs, c(0.05, 0.95)), col = "#AA00CC")
quantile(KLdivObs, 0.992)

# Running the simulation with KL does not give a fundamentally
# different behaviour - since we are just randomly sampling. But KL would be
# more sensitive in case there is biological selection, where the sampling is no
# longer random. If I run the same simulation, with nObs <- 80 but calculating
# KLdiv instead of information, I get a 5% quantile at 0.15 - but the C/L
# frequency swap gives me a KL divergence of 0.18 - this is significant at p =
# 0.008 - (remember, Information is 0 in this case). So that's actually quite a
# nice addition to the toolbox.


# [END]
