# tocID <- "FND-STA-Probability_distribution.R"
#
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the FND-STA-Probability_distribution unit.
#
# Version:  1.6
#
# Date:     2017-10  -  2020-11
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.6    Revise use of data argument to nls()/nlrob()
#           1.5    Extensive additions on Poisson and binomial distributions
#                    and regression fits
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
#TOC>   Section  Title                                                      Line
#TOC> --------------------------------------------------------------------------
#TOC>   1        Introduction                                                 68
#TOC>   2        Three fundamental distributions                             131
#TOC>   2.1        The Poisson Distribution                                  134
#TOC>   2.2        The hypergeometric distribution                           228
#TOC>   2.2.1          Digression: qqplot() and residuals                    376
#TOC>   2.2.2          Fitting functions                                     436
#TOC>   2.2.2.1        Fit a normal distribution (using nls() )              506
#TOC>   2.2.2.2        Fit a normal distribution (using nlrob()) )           525
#TOC>   2.2.2.3        Extreme Value distributions: Gumbel                   552
#TOC>   2.2.2.4        Extreme Value distributions: Weibull                  579
#TOC>   2.2.2.5        Logistic distribution                                 621
#TOC>   2.2.2.6        Log-Logistic distribution                             650
#TOC>   2.2.2.7        Fitting a negative binomial distribution              679
#TOC>   2.2.2.8        Fitting a binomial distribution                       732
#TOC>   2.3        The uniform distribution                                  844
#TOC>   2.4        The Normal Distribution                                   864
#TOC>   3        quantile-quantile comparison                                905
#TOC>   3.1        qqnorm()                                                  915
#TOC>   3.2        qqplot()                                                  981
#TOC>   4        Quantifying the difference                                  998
#TOC>   4.1        Chi2 test for discrete distributions                     1032
#TOC>   4.2        Kullback-Leibler divergence                              1124
#TOC>   4.2.1          An example from tossing dice                         1135
#TOC>   4.2.2          An example from lognormal distributions              1257
#TOC>   4.3        Continuous distributions: Kolmogorov-Smirnov test        1300
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
# 256/6091 ~= 4.2%. How many transcription factors do we expect in a sample of
# 250 genes - like, for example, the top 250 differentially expressed genes in
# an assay? This is a frequently encountered type of question, converging on the
# information contained in a "list of genes". Once an assay has yielded a list
# of genes, what can we learn from looking at staistical features of its
# membership. In the transcription factor model, the question is: how many
# transcription factors do we expect as members of our list - and did we
# actually observe more or less than that?
#
# It would seem that the the probability of encountering _one_ transcription
# factor among our genes is 256/6091 and therefore the number of transcription
# factors we expect to choose at random is (250 * 256)/6091 i.e. 10.50731 in the
# average over many trials. Let's repeat our experiment a million times and see what we get:

N <- 1000000             # One million trials
genes <- numeric(6091)   # All genes are 0
genes[1:256] <- 1        # TFs  are 1

hAssays <- numeric(N)    # initialize vector
set.seed(112358)
for (i in 1:N) {
  pBar(i, N)
  hAssays[i] <- sum(sample(genes, 250)) # sum of TFs in our sample in this trial
}
set.seed(NULL)

# And the average is:
mean(hAssays)  # 10.50293

# ... which is superbly close to our expectation of 10.50731

# All good - but we won't get 10.5 transcription factors in our assay. We'll
# observe five. Or thirteen. Or thirtysix. Or none at all ... and then we ask
# ourselves: is the number of observed transcription factors significantly
# different from what we would have expected if our experiment identified a
# transcription factor just as likely as it identified any other gene? To answer
# this, we need to consider the probability distribution of possible outcomes of
# our assay. Back to the Poisson distribution. In R it is implemented as dpois()
# and its parameters are: the number of observed events, and the probability of
# observing an event:

dpois(0,  (250 * 256) / 6091)       # Probability of seeing no TFs
dpois(1,  (250 * 256) / 6091)       # Probability of seing one ...
dpois(2,  (250 * 256) / 6091)       # Probability of seing two ...
dpois(3:10, (250 * 256) / 6091)     # Probability of seing from three to ten ...
sum(dpois(0:4, (250 * 256) / 6091)) # Probability of seeing four or less ...

sum(dpois(0:250, (250*256)/6091))   # The sum over all possibilities is (for
                                    # any probability distribution) exactly one.

# Lets plot these probabilities ...
nMax <- 28

x <- dpois(0:nMax, (250 * 256) / 6091)

barMids <- barplot(x,
                   col = "#FCF3CF",
                   names.arg = as.character(0:nMax),
                   cex.names = 0.5,
                   axes = TRUE,
                   ylim = c(0, 0.15),
                   xlab = "# of TF in set",
                   ylab = "p")

# ... and add our simulated assays:

(th <- table(factor(hAssays, levels = 0:nMax))/N)

points(barMids - 0.55, th, type = "s", col = "firebrick")
abline(v = (((250 * 256) / 6091) * (barMids[2] - barMids[1])) + barMids[1],
       col = "#5588FF")  # scale to the correct position in the barplot
legend("topright",
       legend = c("Poisson", "simulated", "expectation"),
       pch = c(22, NA, NA),           # one box, two lines only
       pt.cex = 1.4,                  # larger, to match bar-width better
       pt.bg = c("#FCF3CF", NA, NA),  # bg color only for the box
       lty = c(0, 1, 1),              # no line for the box
       col = c("#000000", "firebrick", "#5588FF"),
       bty = "n")                     # no frame around the legend

# NOTE: The simulation shows us that our expectations about the number of
#       transcription factors in a sample are almost, but not exactly Poisson
#       distributed. Can you figure out where we went wrong? Maybe the
#       Poisson distribution is not quite the correct distribution to
#       model our expectations?


# ==   2.2  The hypergeometric distribution  ===================================

# With that suspicion in mind, we reconsider what we were trying to achieve at
# each point of the Poisson distribution. Assume we have observed seven
# transcription factors in our set of 250 genes. So we asked: what is the
# probability of observing seven? And concluded:

dpois(7,  (250 * 256) / 6091)       # Probability of seeing seven ...

# ... which is the probability of seven events in 250 choices if the underlying
# probability is equal to fraction of transcription factors among genes. But
# wait: we weren't careful in our simulation. Assume that our first observed
# gene was a transcription factor. Is the probability of the next sample the
# same? Not quite: one transcription factor has been observed, 249 more samples
# need to be considered, and there are now 6090 genes to choose from. I.e. the
# mean probability changes with every observation:

(250 * 256) / 6091                 # First sample, a TF: p = 10.50731
(249 * (256 - 1)) / (6091 - 1)     # Second sample:      p = 10.42611 (-0.992 %)

# This is actually noticeable: the mean probability for transcription factors
# drops by about one percent after a transcription factor is observed. But what
# if we would have observed the first transcription factor as the tenth gene?

(239 * (256 - 1)) / (6091 - 10)    # Eleventh sample:    p = 10.0222  (-0.954 %)

# This is getting complicated and we need a different way to think about this if
# we don't want to enumerate all possibilities by hand. (There are far too
# many!) Generally, the probability of observing a certain number of events in a
# series is the number of ways to realize the desired outcome, divided by the
# number of possible outcomes. So, if we code transcription factors by "1" and
# other genes by "0", seven transcription factors could be

#   01010100100100100100000000000000000000000 ..., or
#   01110111100000000000000000000000000000000 ..., or
#   11101001000010000000100000000000000000000 ..., or any other combination.

# But crucially, our number of trials is limited and every "success" changes the
# probability of future successes. This is sampling without replacement! And
# sampling without replacement is modeled by the so-called "hypergeometric
# distribution". This is the big difference: when we are sampling WITH
# replacement, we can model the process with a Poisson distribution. When we are
# sampling WITHOUT replacement, we use a hypergeometric distribution instead.

# Let's first re-run our simulated assays. We put the previous run into the
# vector "hAssays" which I prefixed "h" as in "hypergeometric" because I knew
# what was coming, and that I had sampled WITHOUT replacement because that is
# the default for the sample() function. Accordingly, we call the new samples
# "pAssays", where "p" stands for Poisson:

N <- 1000000             # One million trials
genes <- numeric(6091)   # All genes are 0
genes[1:256] <- 1        # TFs  are 1

pAssays <- numeric(N)    # initialize vector
set.seed(112358)
for (i in 1:N) {
  pBar(i, N)
  pAssays[i] <- sum(sample(genes, 250, replace = TRUE))
}
set.seed(NULL)

# Now the average is:
mean(pAssays)  # 10.50312  which is essentially the same as 10.50293

# And the plot ...

nMax <- 28
barMids <- barplot(dpois(0:nMax, (250 * 256) / 6091),
                   names.arg = as.character(0:nMax),
                   cex.names = 0.5,
                   axes = TRUE,
                   ylim = c(0, 0.15),
                   xlab = "# of TF in set",
                   ylab = "p",
                   col = "#FCF3CF")
abline(v = (((250 * 256) / 6091) * (barMids[2] - barMids[1])) + barMids[1],
       col = "#5588FF")  # scale to the correct position in the barplot
th <- table(factor(hAssays, levels = 0:nMax))/N
points(barMids - 0.55, th, type = "s", col = "firebrick")

# Here are the new assays:
tp <- table(factor(pAssays, levels = 0:nMax))/N
points(barMids - 0.55, tp, type = "s", col = "seagreen")

legend("topright",
       legend = c("Poisson",
                  "no replacement",
                  "with replacement",
                  "expectation"),
       pch = c(22, NA, NA, NA),
       pt.cex = 1.4,
       pt.bg = c("#FCF3CF", NA, NA, NA),
       lty = c(0, 1, 1, 1),
       col = c("#000000", "firebrick", "seagreen", "#5588FF"),
       bty = "n")

# Clearly .. the "correct" simulation, the simulation that is actually
# appropriate for the Poisson distribution, matches the theoretical values
# better. Now let's see how well the hypergeometric distribution matches our
# original simulated assays.
# The function dhyper(x, m, n, k) expects the following parameters:
#   x: the number of observed positive events for which to
#      compute the probability
#   m: the number of positive events in the population
#   n: the number of negative eventsw in the population
#   k: the number of observations (or trials)

    dhyper(0,    256, 6091 - 256, 250)  # Probability of seeing no TFs
    dhyper(1,    256, 6091 - 256, 250)  # Probability of seing one ...
    dhyper(2,    256, 6091 - 256, 250)  # Probability of seing two ...
    dhyper(3:10, 256, 6091 - 256, 250)  # Probability of three to ten ...
sum(dhyper(0:4,  256, 6091 - 256, 250)) # Probability of seeing four or less ...

sum(dhyper(0:250, 256, 6091-256, 250)) # The sum over all possibilities is (for
                                       # any probability distribution)
                                       # exactly one.

# Lets plot these probabilities like we did above ...
nMax <- 28
x <- dhyper(0:nMax,  256, 6091 - 256, 250)

barMids <- barplot(x, col = "#E6FFF6",
                   names.arg = as.character(0:nMax),
                   cex.names = 0.5,
                   axes = TRUE,
                   ylim = c(0, 0.15),
                   xlab = "# of TF in set",
                   ylab = "p")
abline(v = (mean(hAssays) * (barMids[2] - barMids[1])) + barMids[1],
       col = "#5588FF")  # scale to the correct position in the barplot
points(barMids - 0.55, th, type = "s", col = "firebrick")

legend("topright",
       legend = c("Hypergeometric",
                  "no replacement",
                  "expectation"),
       pch = c(22, NA, NA),
       pt.cex = 1.4,
       pt.bg = c("#E6FFF6", NA, NA),
       lty = c(0, 1, 1, 1),
       col = c("#000000", "firebrick", "#5588FF"),
       bty = "n")

# This!
# This is what a correctly simulated distribution looks like.


# ===   2.2.1  Digression: qqplot() and residuals               

# The indication that something was wrong with our simulation versus the
# theoretical expectation came from the observation that the differences between
# the barplot and theory were not quite random: the values near the mode were
# systematically too high, the values in the tails were systematically too low.
# This is a general principle: when you see a SYSTEMATIC deviation between
# simulation and theory, something is wrong with your understanding. Either
# there is a subtle error in how you set up the simulation, or a subtle
# misunderstanding about the requirements for the particular theory to apply. R
# has a general way to examine such differences: qqplot() plots the deviations
# between theory and observation ordered as ranks, i.e. not affected by absolute
# scale. (See here for an accessible introduction to qqplot() and the
# "quantiles" that it uses:
# https://data.library.virginia.edu/understanding-q-q-plots/)

oPar <- par(mfrow = c(1, 2))
qqplot( dpois(0:nMax, (250 * 256) / 6091),
        th,
        xlab = "Poisson distribution",
        ylab = "sampling with replacement",
        pch = 22, bg = "#FCF3CF")
abline(lm(th ~ dpois(0:nMax, (250 * 256) / 6091)), col = "seagreen")

qqplot( dhyper(0:nMax,  256, 6091 - 256, 250),
        th,
        xlab = "hypergeometric distribution",
        ylab = "sampling with replacement",
        pch = 22, bg = "#E6FFF6")
abline(lm(th ~ dhyper(0:nMax,  256, 6091 - 256, 250)), col = "seagreen")

par(oPar)

# Similar information can be obtained from a residual plot, plotting differences
# between prediction and observation:

oPar <- par(mfrow = c(1, 2))
plot( dpois(0:nMax, (250 * 256) / 6091) - as.numeric(th),
      type = "b",
      ylim = c(-0.006, 0.006),
      xlab = "# observations",
      ylab = "Poisson density - obs.",
      pch = 22,
      bg = "#FCF3CF",
      col = "#000000")
abline(h = 0, col = "seagreen")

plot( dhyper(0:nMax,  256, 6091 - 256, 250) - as.numeric(th),
      type = "b",
      ylim = c(-0.006, 0.006),
      xlab = "# observations",
      ylab = "Hypergeometric density - obs.",
      pch = 22,
      bg = "#E6FFF6",
      col = "#000000")
abline(h = 0, col = "seagreen")

par(oPar)


# ===   2.2.2  Fitting functions                                

# Note that there is a further subtle catch: We did not actually ask about seven
# transcription factors! We asked about the probability of seven _different_
# transcription factors - because, implicit in the assumptions I made about the
# assay (and reasonable for most gene lists), we count duplicates only once!
# Does this actually make a difference? We can model the situation by giving
# each transcription factor a name (or a number), sampling at random, and
# counting how many unique() factors we found in our sample:

N <- 1000000             # One million trials
genes <- numeric(6091)   # All genes are initially 0
genes[1:256] <- 1:256    # TFs get a unique "name"

# example sample:
set.seed(11235)
x <- sample(genes, 250, replace = TRUE)  # Pick 250 random genes.
x[x != 0]                                # Note that 189 appears twice.
length(x[x != 0])                        # We picked 8 transcription factors ...
unique(x)                                # but only 7 are unique (plus zero).
length(unique(x)) - 1                    # Ignore the zero.

# Do this a million times
uAssays <- numeric(N)    # initialize vector
set.seed(112358)
for (i in 1:N) {
  pBar(i, N)
  uAssays[i] <- length(unique(sample(genes, 250, replace = TRUE))) - 1
}
set.seed(NULL)

# plot the poisson distribution (with replacement) as our baseline
nMax <- 28
barMids <- barplot(dpois(0:nMax, (250 * 256) / 6091),
                   col = "#FCF3CF",
                   names.arg = as.character(0:nMax),
                   cex.names = 0.5,
                   axes = TRUE,
                   ylim = c(0, 0.15),
                   xlab = "# of TF in set",
                   ylab = "p")

tu <- table(factor(uAssays, levels = 0:nMax))/N
points(barMids - 0.55, tu, type = "s", col = "#EE22FF")

legend("topright",
       legend = c("Poisson",
                  "unique genes"),
       pch = c(22, NA),
       pt.cex = 1.4,
       pt.bg = c("#FCF3CF", NA),
       lty = c(0, 1),
       col = c("#000000", "#EE22FF"),
       bty = "n")

# Clearly, the distribution does not match our model exactly.

# So what is the "correct" distribution that we could apply in this case? There
# may or may not be one readily available. What we can do instead is to use a
# more general model and fit parameters. This takes us to the domain of
# regression analysis and curve fitting. The general approach is as follows:
#  - Decide on a statistical model;
#  - Express it in parametrized form;
#  - Fit the parameters;
#  - Analyze the coefficients;

#  Insight into the underlying process that generated our data can be obtained
#  by analyzing the fitted parameters, or simply plotting the results. Let's
#  look at examples of fits to the sampled distribution above:

# ====  2.2.2.1  Fit a normal distribution (using nls() )         

x <- 0:28
plot(x, tu, type="s")
fit <- nls(tu ~ ( a / (sig*sqrt(2*pi)) ) * exp( (-1/2)*((x-mu)/sig)^2 ),
           data = data.frame(x = 0:28, tu = tu),
           start = c(a = 1, mu = 10, sig = 3))  # starting values

points(x, predict(fit), col = "#CC00CC55", lwd = 2,
       type = "s")

coef(fit)
#         a         mu        sig
# 0.9990932 10.0717835  3.0588930

sum(resid(fit)^2)
# [1] 0.0001613488


# ====  2.2.2.2  Fit a normal distribution (using nlrob()) )      

# There's a bit of an art to chosing starting parameters correctly and if the
# nls() fit does not converge, more robust methods are called for.

if (! requireNamespace("robustbase", quietly = TRUE)) {
  install.packages("robustbase")
}

x <- 0:28
plot(x, tu, type="s")
fit <- robustbase::nlrob(tu ~ ( a / (sig*sqrt(6.2831853072)) ) *
                           exp( (-1/2)*((x-mu)/sig)^2 ),
                         data = data.frame(x = 0:28,
                                           tu = tu),
                         start = c(a = 1, mu = 10, sig = 3))  # starting values

points(x, predict(fit), col = "#CC00CC55", lwd = 2, type = "s")

coef(fit)
#        a        mu       sig
# 1.002162 10.059189  3.071217

sum(resid(fit)^2)
# [1] 0.0001630868


# ====  2.2.2.3  Extreme Value distributions: Gumbel              

# Many processes that involve "best-of" choices are better modelled with
# so-called extreme-value distributions: here is the Gumbel distribution
# from the evd package.

if (! requireNamespace("evd", quietly = TRUE)) {
  install.packages("evd")
}

x <- 0:28
plot(x, tu, type="s")

fit <- robustbase::nlrob(tu ~ evd::dgumbel(x, loc = L, scale = S),
                         data = data.frame(tu = tu, x = 0:28),
                         start = c(L = 7.3, S = 2.82))

points(x, predict(fit), type = "s", col = "#55DD88")

coef(fit)
#        L        S
# 9.322110 2.818266

sum(resid(fit)^2)
# [1] 0.001027447


# ====  2.2.2.4  Extreme Value distributions: Weibull             

# Weibull distributions are common in reliabilty analysis. I found the
# distribution particularly hard to fit as it is quite sensitive to inital
# parameter estimates. https://en.wikipedia.org/wiki/Weibull_distribution

# NOTE: the parameter TH is > X
idx <- 4:28
myX <- idx
myY <- as.numeric(tu)[idx]
plot(myX, myY, type="s")

dWei <- function(x, th, l, k) {
  a <- k/l
  b <- ((x-th)/l)^(k-1)
  c <- exp(-((x-th)/l)^k)
  y <- a * b * c
  return(y)
}

set.seed(112358)
fit <- robustbase::nlrob(y ~ dWei(x, th = TH, l = L, k = K),
                         data = data.frame(y = myY,
                                           x = myX),
                         lower = c(TH = 3.5,
                                   L = 8,
                                   K = 2.5),
                         upper = c(TH = 3.7,
                                   L = 9,
                                   K = 3),
                         method = "mtl")

points(x, predict(fit), col = "#CC00CC55", lwd = 2, type = "s")

coef(fit)
#      TH        L        K
#3.630807 8.573898 2.795116

sum(resid(fit)^2)
# [1] 7.807073e-05


# ====  2.2.2.5  Logistic distribution                            

# Similar to normal distribution, but with heavier tails
# https://en.wikipedia.org/wiki/Logistic_distribution

myX <- 0:28
myY <- as.numeric(tu)
plot(myX, myY, type="s")

dLogi <- function(x, mu, s) {
  y <- (exp(-(x-mu)/s)) / (s * (1+exp(-(x-mu)/s))^2)
  return(y)
}

fit <- robustbase::nlrob(y ~ dLogi(x, mu = M, s = S),
                         data = data.frame(y = myY,
                                           x = myX),
                         start = c(M = 10, S = 3))

points(x, predict(fit), col = "#CC00CC55", lwd = 2, type = "s")

coef(fit)
# M         S
# 10.088968  1.891654

sum(resid(fit)^2)
# [1] 0.0004030595


# ====  2.2.2.6  Log-Logistic distribution                        

# Special case of logistic, often used in survival analysis
# https://en.wikipedia.org/wiki/Log-logistic_distribution

myX <- 0:28
myY <- as.numeric(tu)
plot(myX, myY, type="s")

dLogLogi <- function(x, a, b) {
  y <- ((b/a)*(x/a)^(b-1)) / (1+((x/a)^b))^2
  return(y)
}

fit <- robustbase::nlrob(y ~ dLogLogi(x, a = A, b = B),
                         data = data.frame(y = myY,
                                           x = myX),
                         start = c(A = 9.5, B = 4.8))

points(x, predict(fit), col = "#CC00CC55", lwd = 2, type = "s")

coef(fit)
#         A         B
# 10.310181  5.288385

sum(resid(fit)^2)
# [1] 0.0006343927


# ====  2.2.2.7  Fitting a negative binomial distribution         

# The negative binomial distribution is related to the Poisson distribution.
# Assume you are observing events and and counting how many successes of a
# Bernoulli process (essentially a coin-flip) we encounter before we encounter a
# failure. Unlike the Poisson, it models mean and variance separately and is therefore especially useful for overdispersed functions. (cf. https://en.wikipedia.org/wiki/Negative_binomial_distribution)
#

x <- 0:28
plot(x, tu, type="s")

# Negative binomial
dNB <- function(x, r, p) {
  # Note: r is an integer!
  y <- choose((x + r - 1), (r - 1)) * (1-p)^x * p^r
  return(y)
}

set.seed(112358)
RR <- 104
fit <- robustbase::nlrob(tu ~ dNB(x, r = RR, p = P),
                         data = data.frame(x = x, RR = RR),
                         lower = c(P = 0.01),
                         upper = c(P = 0.99),
                         method = "mtl")

points(x, predict(fit), col = "#CC00CC55", lwd = 2, type = "s")

coef(fit)
#         P
# 0.9100729

sum(resid(fit)^2)
# [1] 0.0005669086

# Nb. the parameter R is not continuous: to optimize it as an integer, we try
# reasonable choices and record the best fit.
N <- 150
R <- numeric(N)
for (thisR in 1:N) {
  set.seed(112358)
  fit <- robustbase::nlrob(y ~ dNB(x, r = thisR, p = P),
                           data = data.frame(x = x,
                                             thisR = thisR),
                           lower = c(P = 0.01),
                           upper = c(P = 0.99),
                           method = "mtl")
  R[thisR] <- sum(resid(fit)^2)
}
plot(R)
which(R == min(R))


# ====  2.2.2.8  Fitting a binomial distribution                  
# The workhorse distribution for Bernoulli proceeses
# cf. https://en.wikipedia.org/wiki/Binomial_distribution


myX <- 0:28
myY <- as.numeric(tu)
plot(myX, myY, type="s")

dBinom <- function(x, s, p) {
  y <- choose(s, x) * (p^x) * ((1-p)^(s-x))
  return(y)
}

fit <- robustbase::nlrob(y ~ dBinom(x, s = S, p = P),
                         data = data.frame(y = myY,
                                           x = myX),
                         start = c(S = 240, P = 256/6091))

points(x, predict(fit), col = "#CC00CC55", lwd = 2, type = "s")

coef(fit)
#            S            P
# 122.77746436   0.08376922

sum(resid(fit)^2)
# [1] 1.114272e-06

# Here we go: near perfect fit. But note the parameters! Can you identify the
# relationship of S and P to our model of choosing unique transcription factors
# in a random sampling process?

# Remember what we started out with:
points(x, dBinom(x, 240, 256/6091),
       col = "#00AAFF55", lwd = 2, type = "s")

# The model has improved a lot from our starting-value guess. But how does an
# effective sample size of 122.78 relate to the actual value of 250 samples? And
# how does a probability for choosing a transcription factor of 0.0838 relate to
# having 256 transcription factors among 6091 genes (a fraction of 0.0420)? Well
# - if you stare at the numbers for a bit you might notice that the fitted
# sample size is about half of what we had, while the fitted probability is
# about twice that. Let's correct for the factor of two, and look at the fit ...

points(x, dBinom(x, coef(fit)["S"] * 2 , coef(fit)["P"] / 2),
       col = "#00FF8899", lwd = 2, type = "s")

# interesting ... that's not so bad ... Essentially what this tells us is: our
# model is very similar to one in which our effective sample size would have
# been about 246 genes, not 250, and the number of our genes would have been
# close to 0.04188461 * 6091 ~= 255 transcription factors, not 256. Close, to -
# but not almost perfect.

# Incidentally: would you have thought that we can
# detect the number of transcription factors in our sample +- 1 if we just
# sample frequently enough?

# Let's examine the effect that applying various factors to size and probability
# have on our distribution:


z <- c(1/8:2 , 1, 2:8, 10, 20, 50)
myCols <- colorRampPalette(c("#00000022",
                             "#00000088",
                             "#DD0000DD",
                             "#0088FF88",
                             "#0088FF55",
                             "#0088FF22"),
                           alpha=TRUE)(length(z))

myLwd <- rep(1.25, length(z))
myLwd[z == 1] <- 2.5 # make the actual fitted value stand out

plot(x, y, type="s", ylim =c(0, 0.22), lwd = 5, col="#FF990077")

for (i in seq_along(z)) {
  points(x, dBinom(x, coef(fit)["S"] * z[i] , coef(fit)["P"] / z[i]),
         col = myCols[i], type = "s")
}

legend("topright",
       legend = sprintf("%3.2f", z),
       cex = 0.5,
       col = myCols,
       lwd = myLwd,
       bty = "n")

# We see how factors less than 1 applied to our actual sample size and fraction
# of transcription factors increasingly narrow the distribution; factors larger
# than one increasingly broaden the distribution but in a much slower manner. A
# factor of one (i.e. half size, double probability) hits the sweet spot and
# recaptures the essence of our observation. Clearly, this is a case of
# fortuitously cancelling errors.


# What is the take-home message here?
#  - The standard Poisson and hypergeometric distributions apply to
#      very specific models of choice processes (with/without replacement);
#  - It is not trivial to choose the right statistical model for any particular
#      real-world specification of a sampling process. Our model of unique
#      choices is neither Poisson nor hypergeometric distributed. Once we have
#      identified the parameters of a binomial distribution that models the
#      process perfectly, we nevertheless note that we don't seem to be
#      able to interpret the parameters easily.
#  - Stochastic modeling allows us to validate whether the model is correct,
#      but in case there are discrepancies it is not obvious what might be
#      a more appropriate model and how to parametrize it.
#  - If building an explicit stochastic model is not possible, we'll have to
#      to do the best we can, in this case that would mean: choose a more
#      general distribution and parametrize it.


# ==   2.3  The uniform distribution  ==========================================

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


# ==   2.4  The Normal Distribution  ===========================================

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
# (https://stats.stackexchange.com/questions/113692)
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


# ==   4.3  Continuous distributions: Kolmogorov-Smirnov test  =================

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
