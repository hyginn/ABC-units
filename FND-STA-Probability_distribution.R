# tocID <- "FND-STA-Probability_distribution.R"
#
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the FND-STA-Probability_distribution unit.
#
# Version:  1.4
#
# Date:     2017-10  -  2020-09
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.4    2020 Maintenance
#           1.3    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout,
#           1.2    Update set.seed() usage
#           1.1    Corrected empirical p-value
#           1.0    First code live version
#
# TODO:
#           Add tasks
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
#TOC>   Section  Title                                                         Line
#TOC> -----------------------------------------------------------------------------
#TOC>   1        Introduction                                                    54
#TOC>   2        Three fundamental distributions                                117
#TOC>   2.1        The Poisson Distribution                                     120
#TOC>   2.2        The uniform distribution                                     174
#TOC>   2.3        The Normal Distribution                                      194
#TOC>   3        quantile-quantile comparison                                   235
#TOC>   3.1        qqnorm()                                                     245
#TOC>   3.2        qqplot()                                                     311
#TOC>   4        Quantifying the difference                                     328
#TOC>   4.1        Chi2 test for discrete distributions                         363
#TOC>   4.2        Kullback-Leibler divergence                                  454
#TOC>   4.2.1          An example from tossing dice                             465
#TOC>   4.2.2          An example from lognormal distributions                  588
#TOC>   4.3        Kolmogorov-Smirnov test for continuous distributions         631
#TOC> 
#TOC> ==========================================================================


# =    1  Introduction  ========================================================

# The space of possible outcomes of events is called a probability distribution
# and the properties of probability distributions are crucial to our work. Many
# distributions like the (discrete) Poisson, or the (continuous) Normal
# distribution have been characterized in great detail, their behaviour is well
# understood, and they are useful for countless applications - but we also need
# to able to work with ad hoc distributions that have never been seen before.

# Let's get a few facts about probability distributions out of the way:

# The "support" of a probability distribution is the range of outcomes that have
# a non-zero probability. The "domain" of a probability distribution is the
# range of probabilities that the distribution can take over its support. Think
# of this as the ranges on the x- and y-axis respectively. Thus the distribution
# can be written as p = f(x).

# The integral over a probability distribution is always 1. This means: the
# distribution reflects the situation that an event does occur, any event, but
# there is not "no event".

# R's inbuilt probability functions always come in four flavours:
# d...  for "density": this is the probability density function (p.d.f.),
#         the value of f(x) at x.
# p...  for "probability": this is the cumulative distribution function
#         (c.d.f.). It is 0 at the left edge of the support, and 1 at
#         the right edge.
# q... for "quantile": The quantile function returns the x value at which p...
#         takes a requested value.
# r... for "random": produces random numbers that are distributed according
#         to the p.d.f.

# To illustrate with the "Normal Distribution" (Gaussian distribution):

# 1000 normally distributed values with default parameters: mean 0, sd 1.
r <- rnorm(1000)

# pastel green: histogram of 1000 random samples
hist(r,
     freq = FALSE,
     breaks = 30,
     xlim = c(-4, 4),
     ylim = c(0, 1),
     main = "Normal Distribution",
     xlab = "x",
     ylab = "f(x)",
     col = "#E6FFF6")

# 100 equally spaced point along x
x <- seq(-4, 4, length.out = 100)

# black: c. d. f. along x. Note that this asymptotically reaches 1
points(x, pnorm(x), type = "l")

# dark red: p. d. f. along x
points(x, dnorm(x), type = "l", lwd = 2, col="firebrick")

# purple: 1% and 99% quantiles
abline(v = qnorm(c(0.01, 0.99)), lwd = 0.5, col = "#CCAAFF")

# Study this plot well and familiarize yourself with the terms.


# =    2  Three fundamental distributions  =====================================


# ==   2.1  The Poisson Distribution  ==========================================

# The Poisson distribution is a discrete probability distribution that
# characterizes how many events we expect among a given number of observations
# if the mean probability of an event is known. Assume we know that there are
# 256 transcription factors among the 6091 protein-coding genes of yeast, then
# the probability of picking a transcription factor at random from all ORFs is
# 256/6091 ~= 4.2%. How many do we expect if we look e.g. at 250 differentially
# expressed genes? This means the mean number of transcription factors  we would
# expect in that sample of differentially expressed genes is (250 * 256)/6091.

dpois(0,  (250 * 256) / 6091)       # Probability of seeing no TFs
dpois(1,  (250 * 256) / 6091)       # Probability of seing one ...
dpois(2,  (250 * 256) / 6091)       # Probability of seing two ...
dpois(3:10, (250 * 256) / 6091)     # Probability of seing from three to ten ...
sum(dpois(0:4, (250 * 256) / 6091)) # Probability of seeing four or less ...

# Lets plot this
N <- 25
x <- dpois(0:N, (250 * 256) / 6091)
names(x) <- as.character(0:N)
midPoints <- barplot(x, col = "#E6FFF6",
             axes = TRUE,
             ylim = c(0, 0.15),
             xlab = "# of TF in set",
             ylab = "p")

# Confirm that our understanding of dpois() is correct, by simulating actual
# trials:

N <- 1000
genes <- numeric(6091)   # All genes are 0
genes[1:256] <- 1        # TFs  are 1

x <- numeric(N)          # initialize vector
set.seed(112358)
for (i in 1:N) {
  x[i] <- sum(sample(genes, 250)) # sum of TFs in our sample in this trial
}
set.seed(NULL)

(t <- table(x)/N)

# Add these values to the plot
y <- numeric(26)                     # initialize vector with 26 slots
y[as.numeric(names(t)) + 1] <- t     # put the tabled values there (index + 1)
points(midPoints - 0.55, y, type = "s", col = "firebrick")
legend("topright",
       legend = c("theoretical", "simulated"),
       pch = c(22, 22),
       pt.bg = c("#E6FFF6", "firebrick"),
       bty = "n")


# ==   2.2  The uniform distribution  ==========================================

# The uniform distribution has the same probability over its entire support. R's
# runif() function takes the desired number, the min and the max as arguments.
# Let's plot a histogram of a million values between -1 and 1 to demonstrate.

hist(runif(1e6, -1, 1), breaks = 20, col = "#FFE6F6")
abline(h = 1e6/20, col="firebrick")

# One of the important uses of the uniform distribution is to write conditional
# expressions that are TRUE with a given probability. For example, to get TRUE
# values with a probability of 1/3, we pick a random number between 0 and 1, and
# ask whether the number is smaller than 1/3. Example:

runif(10) < 1/3

#confirm with a million trials
sum(runif(1e6) < 1/3)/1e6  # should be close to 0.33333...


# ==   2.3  The Normal Distribution  ===========================================

# The king of probability distributions. Why? That's because of the Central
# Limit Theorem (CLT) that essentially says that a process that is subject to
# small fluctuations will tend to normally distributed outcomes ... and in
# biology literally everything is subject to small fluctuations.

# Lets simulate: Lets create a vector of 12345 numbers, choose samples of size
# 77 from that vector and calculate the mean. We do that 9999 times. I am just
# using odd numbers so it is clear in the code which number represents what.

x <- runif(12345)
v <- numeric(9999)
for (i in 1:length(v)) {
  v[i] <- mean(sample(x, 77))
}

hist(v, breaks = 20, col = "#F8DDFF")

# Let's try this again, this time sampling from an exponential distribution ...
x <- rexp(12345)
v <- numeric(9999)
for (i in 1:length(v)) {
  v[i] <- mean(sample(x, 77))
}

hist(v, breaks = 20, col = "#F8DDFF")

# ... or from a t-distribution ...
x <- rt(12345, df = 3)
v <- numeric(9999)
for (i in 1:length(v)) {
  v[i] <- mean(sample(x, 77))
}

hist(v, breaks = 20, col = "#F8DDFF", freq = FALSE)

# The outcomes all give normal distributions, regardless what the details of our
# original distribution were!


# =    3  quantile-quantile comparison  ========================================


# Once we have observed a distribution of outcomes, the next task is often to
# quantify whether the values correspond to a known distribution. This can be
# done with quantile-quantile plots that plot the quantile values of one
# distribution against those of the other. Identically distributed quantiles lie
# on a straight line.


# ==   3.1  qqnorm()  ==========================================================

# The functions qqnorm() and qqline() perform this
# comparison with the normal distribution.

set.seed(112358)
x <- rnorm(100, mean=0, sd=1) # 100 normally distributed values
set.seed(NULL)

qqnorm(x)
qqline(x, col = "seagreen")

# Our variables in x appear normally distribute - which is not surprising since
# we produced them with rnorm(). What about the kind of sampling we did above?
# Let's test whether samples from an exponential distribution are normally
# distributed (previously we just visually inspected the histogram.)

# Create a vector of sample means from the exponential distribution; use
# only a few samples for the mean
x <- rexp(12345)
v <- numeric(999)

set.seed(112358)
for (i in 1:length(v)) {
  v[i] <- mean(sample(x, 12))
}
set.seed(NULL)

qqnorm(v)
qqline(v, col = "turquoise") # normal

#try with more samples
for (i in 1:length(v)) {
  v[i] <- mean(sample(x, 77))
}
qqnorm(v)
qqline(v, col = "steelblue") # normal

#try with many samples
for (i in 1:length(v)) {
  v[i] <- mean(sample(x, 374))
}
qqnorm(v)
qqline(v, col = "plum") # normal

# Exactly as the CLT predicts, the more often we sample - we tried
# 12, 77, and 374 samples - the more "normal" the distribution looks.

# What does a distribution look like that is NOT normal?
# Let's try simulating an "Extreme value distribution". This type of
# distribution appears if we choose the max() of a sample

set.seed(112358)
rEVD <- numeric(9999)
for (i in seq_along(rEVD)) {
  rEVD[i] <- max(rnorm(100))
}
set.seed(NULL)

hist(rEVD, breaks = 20, col = "orchid")
# Note the long tail on the right!

qqnorm(rEVD)
qqline(rEVD, col = "orchid") # Definitely not "normal"!


# ==   3.2  qqplot()  ==========================================================

# qqplot() works like qqnorm(), except that we compare two arbitrary
# distribtutions, rather than one distribution against the normal distribution
# as we did before. Sample against sample.

x <- seq(0, 4, length.out = 100)

dl <- dlnorm(x) # log-normal distribution
dg <- dgamma(x, shape=1.5)   # gamma distribution

plot(dl, type="l")
points(dg, type="l", col = "maroon")

qqplot(dl, dg) # clearly not equal


# =    4  Quantifying the difference  ==========================================


# Quantile-quantile plots give us a visual estimate, but how do we quantify the
# difference between distributions? Let's compare two types of extreme-value
# distributions, a lognormal distribution with the same distribution that is
# slightly shifted, and three gamma distributions with three different shape
# parameters. Let's define and visualize the distributions first, to see what
# we are comparing.

x <- seq(0, 4, length.out = 100)

set.seed(112358)
dl1 <- dlnorm(x) # log-normal distribution
dl2 <- dlnorm(x - 0.25) # log-normal distribution, shifted right (a bit)
dg1.2 <- dgamma(x, shape=1.2)   # three gamma distributions with...
dg1.5 <- dgamma(x, shape=1.5)   # ...wider, and wider...
dg1.9 <- dgamma(x, shape=1.9)   # ...peak
set.seed(NULL)

myCols <- c("black", "grey", "maroon", "turquoise", "steelblue")

plot(dl1, type="l", lwd=2)   # visualize the distributions
points(dl2, type="l", col = myCols[2])
points(dg1.2, type="l", col = myCols[3])
points(dg1.5, type="l", col = myCols[4])
points(dg1.9, type="l", col = myCols[5])
legend("topright",
       legend = c("dl1", "dl2", "dg1.2", "dg1.5", "dg1.9"),
       lty = 1,
       lwd = c(2, 1, 1, 1, 1),
       col = myCols,
       bty = "n")


# ==   4.1  Chi2 test for discrete distributions  ==============================

# The chi2 test can be used to compare discrete distributions - even though
# it is not ideal for this purpose.
# (https://stats.stackexchange.com/questions/113692/test-identicality-of-discrete-distributions)
#
# Let's draw 100 random variables according to our target distributions
N <- 100
set.seed(112358)
rL1   <- rlnorm(N) # log-normal distribution
rL2   <- rlnorm(N, meanlog = 0.25) # log-normal distribution, shifted right
rG1.2 <- rgamma(N, shape=1.2)   # three gamma distributions with...
rG1.5 <- rgamma(N, shape=1.5)   # ...wider, and wider...
rG1.9 <- rgamma(N, shape=1.9)   # ...peak
set.seed(NULL)

maxX <- max(c(rL1, rL2, rG1.2, rG1.5, rG1.9))

myBreaks <- seq(0, 5, length.out = 10) # 9 intervals from 0 to 5...
myBreaks <- c(myBreaks, maxX) # ... and one that contains the outliers

# It's easy to plot a histogram of one set of random deviates...
hist(rG1.5, breaks = myBreaks, col = myCols[4])

# ... but basic R has no inbuilt function to stack histogram bars side-by-side.
# We use the multhist() function in the plotrix package: check out the
# package information - plotrix has _many_ useful utilities to enhance
# plots or produce informative visualizations.

if (! requireNamespace("plotrix", quietly = TRUE)) {
  install.packages("plotrix")
}
# Package information:
#  library(help = plotrix)       # basic information
#  browseVignettes("plotrix")    # available vignettes
#  data(package = "plotrix")     # available datasets


h <- plotrix::multhist(list(rL1, rL2, rG1.2, rG1.5, rG1.9 ),
                       breaks = myBreaks,
                       col = myCols)
legend("topright",
       legend = c("rL1", "rL2", "rG1.2", "rG1.5", "rG1.9"),
       pch=15,
       col = myCols,
       bty="n")

# We have assigned the output of multihist to the variable h, which now
# contains the densities, midpoints, breaks etc. of the histogram ...
# as well as the matrix of counts:
h[[2]]

# But using hist() or multhist() for binning is a bit of a hack - even though
# it works, of course. The "real" R function to bin data is cut(). cut()
# sorts the values of a vector into bins, and returns the bin-labels as
# integers.

x <- cut(rL1, breaks = myBreaks, labels = FALSE)
table(x)

countsL1   <- table(cut(rL1  , breaks = myBreaks, labels = FALSE))
countsL2   <- table(cut(rL2  , breaks = myBreaks, labels = FALSE))
countsG1.2 <- table(cut(rG1.2, breaks = myBreaks, labels = FALSE))
countsG1.5 <- table(cut(rG1.5, breaks = myBreaks, labels = FALSE))
countsG1.9 <- table(cut(rG1.9, breaks = myBreaks, labels = FALSE))


chisq.test(countsL1, countsL2)

# Note that this use of the chi2 test does not have very much power, since it
# does not assume the values for the bins are ordered, but takes them to
# be independent. Also, chisq.test() complains about the Chi-squared
# approximation to be possibly incorrect because some of the counts are small.
# In this case, rather than rely on the in-built chi-square table, we can do
# an explicit simulation to estimate the p-value for the null hypothesis.
#
chisq.test(countsL1, countsL2, simulate.p.value = TRUE, B = 10000)

# Note that the probability that the samples came from the same distribution is
# quite high. We don't seem to be able to distinguish l1 and l2 with the chi2
# test.

# Let's calculate this for the other distribution too:

chisq.test(countsL1, countsG1.2, simulate.p.value = TRUE, B = 10000)
chisq.test(countsL1, countsG1.5, simulate.p.value = TRUE, B = 10000)
chisq.test(countsL1, countsG1.9, simulate.p.value = TRUE, B = 10000)

# As a result, we would conclude that none of these distributions are
# significantly different.

# ==   4.2  Kullback-Leibler divergence  =======================================

# For discrete probability distributions, there is a much better statistic, the
# Kullback-Leibler divergence (or relative entropy). It is based in information
# theory, and evaluates how different the matched pairs of outcome categories
# are. Its inputs are the probability mass functions (p.m.f.) of the two
# functions to be compared. A probability mass function is the probability of
# every outcome the process can have. Kullback-Leibler divergence therefore can
# be applied to discrete distributions. But we need to talk a bit about
# converting counts to p.m.f.'s.

# ===   4.2.1  An example from tossing dice                        

#  The p.m.f of an honest die is (1:1/6, 2:1/6, 3:1/6, 4:1/6, 5:1/6, 6:1/6). But
#  there is an issue when we convert sampled counts to frequencies, and estimate
#  probabilities with these frequencies. The problem is: by chance we might not
#  have observed a particular outcome! Consider this example:

set.seed(47)
N <- 20
(counts <- table(sample(1:6, N, replace = TRUE)))
set.seed(NULL)

# We have not observed a "2"!
#
# Now, if we convert the counts to frequencies, and assume these to be
# probabilities, we get an observed p.m.f:
pmf <- numeric(6)
for (i in seq_along(counts)) {
  outcome <- as.integer(names(counts)[i])
  pmf[outcome] <- counts[i]/N
}
names(pmf) <- 1:6
pmf

# This means we assign a probability of 0 to the outcome "2" - because we
# haven't observed it. But how do we compare it then? What does it mean when we
# should have 1/6 of "2"s but we have none whatsoever? Taken at face value it
# means we are comparing apples to oranges - essentially a five-faced die with a
# six-faced die - and therefore divergence is infinitely large. Obviously, that
# was just a problem caused by our limited sampling, but the question remains:
# how do we account for missing data? There are several solutions, for example,
# for ordered data one could substitute the average values of the two bracketing
# outcomes. But a simple and quite robust solution is to add "pseudocounts".
# This is called adding a Laplace prior, or a Jeffreys prior: in our case,
# simply add 0.5 to every category.

# pmf of an honest die
pmfHD <- rep(1/6, 6)
names(pmfHD) <- 1:6
pmfHD

# pmf, estimated from our sampled counts with and without pseudocounts

pmfPC <- numeric(6) + 0.5
for (i in seq_along(counts)) {
  outcome <- as.integer(names(counts)[i])
  pmfPC[outcome] <- counts[i] + pmfPC[outcome]
}
pmfPC <- pmfPC / sum(pmfPC)
names(pmfPC) <- 1:6
pmf     # before
pmfPC   # after

# The definition of a function that takes samples and returns a pmf
# with pseudocounts is straightforward, but requires a bit of juggling
# the observed values into the right slots. Here's what this can look like:

pmfPC <- function(cts, nam) {
  # Convert counts to a p.m.f, adding a pseudocount of 0.5 to all categories
  # of outcomes (Jeffreys prior).
  # Parameters:
  #    cts  num   raw counts
  #    nam        names of the categories, converted to char
  # Value   a probability mass function

  nam <- as.character(nam)
  if (!all(names(cts) %in% nam)) {
    stop("PANIC: names in \"cts\" do not match \"nam\"!")
  }
  pmf <- numeric(length(nam))
  names(pmf) <- nam
  pmf[names(cts)] <- cts[names(cts)]
  pmf <- pmf + 0.5    # add pseudocounts
  return(pmf/sum(pmf))
}



# The definition of the Kullback-Leibler divergence itself is quite simple
# actually: for two distributions, p and q it is sum(p * log( p / q))

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

# Now we can calculate KL-divergences for our six-faced die example:

KLdiv(rep(1/6, 6), pmfPC(counts, 1:6)) # p.m.f. of an honest die vs. our
                                       # actual counts

# Is that a lot? Let's look at the distribution of KL-divergences in our
# six-sided die scenario: p.m.f. of 20 samples each versus the theoretical
# p.m.f.

N <- 1000
nSample <- 20
divs <- numeric(N)
for (i in 1:N) {
  x <- table(sample(1:6, nSample, replace = TRUE))
  divs[i] <- KLdiv(rep(1/6, 6), pmfPC(x, 1:6))
}

hist(divs,
     col = "whitesmoke",
     xlab = "Kullback-Leibler divergence",
     main = sprintf("%d samples of honest die", nSample))
abline(v = KLdiv(rep(1/6, 6), pmfPC(counts, 1:6)), col="firebrick")

# We see that our example sample (KL-divergence marked with the red line) is
# somewhat but not drastically atypical.


# ===   4.2.2  An example from lognormal distributions             

# We had compared a set of lognormal and gamma distributions above, now we
# can use KL-divergence to quantify their similarity:

pmfL1 <- pmfPC(countsL1, nam = 1:10)
pmfL2 <- pmfPC(countsL2, nam = 1:10)
KLdiv(pmfL1, pmfL2)  # 0.1087

# To evaluate what this number means, we can run a simple simulation: we create
# random samples according to the rL1 distribution, calculate the Kullback
# Leibler divergence with countsL1, and compare the distribution we get with the
# value we observed as the difference with discL2. Essentially, this tells us
# the probability that countsL2 is actually a sample from the L1 function.
# Here we go:

N <- 1000
divs <- numeric(N)

set.seed(31416)
for (i in 1:N) {
  x <- rlnorm(100) # get random samples from our reference distributions
  y <- table(cut(x, breaks = myBreaks, labels = FALSE))  # tabulate counts
  q <- pmfPC(y, nam = 1:10)  # convert to p.m.f. with pseudocounts
  divs[i] <- KLdiv(pmfL1, q)     # calculate Kullback-Leibler divergence
}
set.seed(NULL)

hist(divs,
     col = "thistle",
     xlab = "Kullback-Leibler divergence",
     main = sprintf("Samples from shifted lnorm()"))
abline(v = KLdiv(pmfL1, pmfL2), col="firebrick")

# How many KL-divergences were less than the difference we observed?
sum(divs < KLdiv(pmfL1, pmfL2)) # 933

# Therefore the empirical p-value that the samples came from the same
# distribution is only 100 * ((N - 933) + 1) / (N + 1) (%) ... 6.8%. You see
# that this gives a much more powerful statistical approach than the chi2 test
# we used above.


# ==   4.3  Kolmogorov-Smirnov test for continuous distributions  ==============

# The Kolmogorov-Smirnov (KS) test is meant for continuous distributions, i.e.
# the probability it calculates assumes that the function values are all
# different. In the case of random samples, that is a reasonable assumption. KS
# is a two-sample goodness of fit test, i.e. it has a null-hypothesis that the
# samples were taken from the same (unknown) distribution, and then asks whether
# this is likely, given the actual values. It is sensitive both to differences
# in location (position of the mean), and differences in shape. R has the
# ks.test() function to calculate it.

ks.test(dl1, dl2)     # vs. shifted log-normal.
ks.test(dl1, dg1.2)   # vs. the three gamma distributions
ks.test(dl1, dg1.5)
ks.test(dl1, dg1.9)
ks.test(dg1.5, dg1.9) # the two last gammas against each other

# (The warnings about the presence of ties comes from the 0's in our function
# values). The p-value that these distributions are samples of the same
# probability distribution gets progressively smaller.



# [END]
