# tocID <- "FND-MAT-Graphs_and_networks.R"
#
# ---------------------------------------------------------------------------- #
#  PATIENCE  ...                                                               #
#    Do not yet work wih this code. Updates in progress. Thank you.            #
#    boris.steipe@utoronto.ca                                                  #
# ---------------------------------------------------------------------------- #
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the FND-MAT-Graphs_and_networks unit.
#
# Version:  1.2
#
# Date:     2017  10  -  2019  01
# Author:   Boris Steipe (boris.steipe@utoronto.ca)
#
# Versions:
#           1.2    Change from require() to requireNamespace(),
#                      use <package>::<function>() idiom throughout
#           1.1    Update set.seed() usage
#           1.0    First final version for learning units.
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
#TOC>   Section  Title                                        Line
#TOC> ------------------------------------------------------------
#TOC>   1        Review                                         50
#TOC>   2        DEGREE DISTRIBUTIONS                          204
#TOC>   2.1        Random graph                                210
#TOC>   2.2        scale-free graph (Barabasi-Albert)          258
#TOC>   2.3        Random geometric graph                      323
#TOC>   3        A CLOSER LOOK AT THE igraph PACKAGE           445
#TOC>   3.1        Basics                                      448
#TOC>   3.2        Components                                  520
#TOC>   4        RANDOM GRAPHS AND GRAPH METRICS               539
#TOC>   4.1        Diameter                                    576
#TOC>   5        GRAPH CLUSTERING                              645
#TOC>
#TOC> ==========================================================================


# =    1  Review  ==============================================================

# This tutorial covers basic concepts of graph theory and analysis in R. Make
# sure you have pulled the latest version of the project from the GitHub
# repository, and that you have typed init() to load some utility functions and
# data.

# Let's explore some of the basic ideas of graph theory by starting with a small
# random graph.


# To begin let's write a little function that will create random "gene" names;
# there's no particular purpose to this other than to make our graphs look a
# little more "biological" ...
makeRandomGenenames <- function(N) {
  nam <- character()
  while (length(nam) < N) {
    a <- paste0(c(sample(LETTERS, 1), sample(letters, 2)),
                collapse="") # one uppercase, two lowercase letters
    n <- sample(1:9, 1)      # one number
    nam[length(nam) + 1] <- paste(a, n, sep="") # store in vector
    nam <- unique(nam)   # delete if this was a duplicate
  }
  return(nam)
}

N <- 20

set.seed(112358)                       # set RNG seed for repeatable randomness
(Nnames <- makeRandomGenenames(N))
set.seed(NULL)                         # reset the RNG

# One way to represent graphs in a computer is as an "adjacency matrix". In this
# matrix, each row and each column represents a node, and the cell at the
# intersection of a row and column contains a value/TRUE if there is an edge,
# 0/FALSE otherwise.

# Let's create an adjacency matrix for random graph: let's say a pair of nodes
# has probability p <- 0.1 to have an edge, and our graph is symmetric , i.e. it
# is an undirected graph, and it has neither self-edges, i.e. loops, nor
# multiple edges between the same nodes, i.e. it is a "simple" graph. We use our
# the Nnames vector as node labels.

makeRandomAM <- function(nam, p = 0.1) {
  # Make a random adjacency matrix for a set of nodes with edge probability p
  # Parameters:
  #   nam: a character vector of unique node names.
  #   p:   probability that a random pair of nodes will have an edge.
  #
  # Value: an adjacency matrix for a simple, undirected graph
  #

  N <- length(nam)
  AM <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
  rownames(AM) <- nam
  colnames(AM) <- nam
  for (iRow in 1:(N-1)) { # Note how we make sure iRow != iCol - this prevents
                          # loops
    for (iCol in (iRow+1):N) {
      if (runif(1) < p) {     # runif() creates uniform random numbers
                              # between 0 and 1. The expression is TRUE with
                              # probability p. if it is TRUE ...
        AM[iRow, iCol] <- 1   # ... record an edge for the pair (iRow, iCol)
      }
    }
  }
  return(AM)
}

set.seed(112358)                       # set RNG seed for repeatable randomness
(myRandAM <- makeRandomAM(Nnames, p = 0.09))
set.seed(NULL)                         # reset the RNG


# Listing the matrix is not very informative - we should plot this graph. The
# standard package for work with graphs in r is "igraph". We'll go into more
# details of the igraph package a bit later, for now we just use it to plot:

if (! requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
# Package information:
#  library(help = igraph)       # basic information
#  browseVignettes("igraph")    # available vignettes
#  data(package = "igraph")     # available datasets


myG <- igraph::graph_from_adjacency_matrix(myRandAM, mode = "undirected")

set.seed(112358)                       # set RNG seed for repeatable randomness
                                       # calculate layout coordinates
myGxy <- igraph::layout_with_graphopt(myG,
                                      charge=0.0012)
set.seed(NULL)                         # reset the RNG


# The igraph package adds its own function to the collection of plot()
# functions; R makes the selection which plot function to use based on the class
# of the object that we request to plot. This plot function has parameters
#  layout - the x,y coordinates of the nodes;
#  vertex.color - which I define to color by node-degree
#  vertex size - which I define to increase with node-degree
#  vertex.label - which I set to combine the names of the vertices of the
#                 graph - names(V(iG)) - with the node degree - degree(iG).
# See ?igraph.plotting for the complete list of parameters


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myG,
     layout = myGxy,
     rescale = FALSE,
     xlim = c(min(myGxy[,1]) * 0.99, max(myGxy[,1]) * 1.01),
     ylim = c(min(myGxy[,2]) * 0.99, max(myGxy[,2]) * 1.01),
     vertex.color=heat.colors(max(igraph::degree(myG)+1))[igraph::degree(myG)+1],
     vertex.size = 1600 + (300 * igraph::degree(myG)),
     vertex.label = sprintf("%s(%i)", names(igraph::V(myG)), igraph::degree(myG)),
     vertex.label.family = "sans",
     vertex.label.cex = 0.7)
par(oPar)  # reset plot window


# The simplest descriptor of a graph are the number of nodes, edges, and the
# degree-distribution. In our example, the number of nodes was given: N; the
# number of edges can easily be calculated from the adjacency matrix. In our
# matrix, we have entered 1 for every edge. Thus we simply sum over the matrix:
sum(myRandAM)

# Is that what you expect?

# What about the degree distribution? We can get that simply by summing over the
# rows and summing over the columns and adding the two vectors.

rowSums(myRandAM) +  colSums(myRandAM) # check this against the plot!

# The function degree() gives the same values
igraph::degree(myG)

# Let's  plot the degree distribution in a histogram:
degG <- igraph::degree(myG)
brk <- seq(min(degG)-0.5, max(degG)+0.5, by=1)  # define histogram breaks
hist(degG, breaks=brk, col="#A5CCF5",
     xlim = c(-1,8), xaxt = "n",
     main = "Node degrees", xlab = "Degree", ylab = "Number")  # plot histogram
axis(side = 1, at = 0:7)

# Note: I don't _have_ to define breaks, the hist() function usually does so
# quite well, automatically. But for this purpose I want the columns of the
# histogram to represent exactly one node-degree difference.

# A degree distribution is actually quite an important descriptor of graphs,
# since it is very sensitive to the generating mechanism. For biological
# networks, that is one of the key questions we are interested in: how was the
# network formed?

# =    2  DEGREE DISTRIBUTIONS  ================================================

# Let's simulate a few graphs that are a bit bigger to get a better sense of
# their degree distributions:
#

# ==   2.1  Random graph  ======================================================

N <- 200

set.seed(31415927)                     # set RNG seed for repeatable randomness
my200AM <- makeRandomAM(as.character(1:N), p = 0.015)
set.seed(NULL)                         # reset the RNG

myG200 <- igraph::graph_from_adjacency_matrix(my200AM, mode = "undirected")
myGxy <- igraph::layout_with_graphopt(myG200, charge=0.0001) # calculate layout
                                                             # coordinates

oPar <- par(mar= rep(0,4))             # Turn margins off, save graphics state
plot(myG200,
     layout = myGxy,
     rescale = FALSE,
     xlim = c(min(myGxy[,1]) * 0.99, max(myGxy[,1]) * 1.01),
     ylim = c(min(myGxy[,2]) * 0.99, max(myGxy[,2]) * 1.01),
     vertex.color=heat.colors(max(igraph::degree(myG200)+1))[igraph::degree(myG200)+1],
     vertex.size = 150 + (60 * igraph::degree(myG200)),
     vertex.label = NA)
par(oPar)                              # restore graphics state

# This graph has thirteen singletons and one large, connected component. Many
# biological graphs look approximately like this.

# Calculate degree distributions
dg <- igraph::degree(myG200)
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5F5CC",
     xlim = c(-1,11), xaxt = "n",
     main = "Node degrees", xlab = "Degree", ylab = "Number")  # plot histogram
axis(side = 1, at = 0:10)


# Note the pronounced peak of this distribution: this is not "scale-free".
# Here is the log-log plot of frequency vs. degree-rank ...

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#A5F5CC",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a random network")

# ... which shows us that this does NOT correspond to the single-slope linear
# relationship that we expect for a "scale-free" graph.

# ==   2.2  scale-free graph (Barabasi-Albert)  ================================

# What does one of those intriguing "scale-free" distributions look like? The
# iGraph package has a function to make random graphs according to the
# Barabasi-Albert model of scale-free graphs. It is: sample_pa(), where pa
# stands for "preferential attachment". Preferential attachment is one type of
# process that will yield scale-free distributions.

N <- 200

set.seed(31415927)                     # set RNG seed for repeatable randomness
GBA <- igraph::sample_pa(N, power = 0.8, directed = FALSE)
set.seed(NULL)                         # reset the RNG

GBAxy <- igraph::layout_with_graphopt(GBA, charge=0.0001)

oPar <- par(mar= rep(0,4))             # Turn margins off, save graphics state
plot(GBA,
     layout = GBAxy,
     rescale = FALSE,
     xlim = c(min(GBAxy[,1]) * 0.99, max(GBAxy[,1]) * 1.01),
     ylim = c(min(GBAxy[,2]) * 0.99, max(GBAxy[,2]) * 1.01),
     vertex.color=heat.colors(max(igraph::degree(GBA)+1))[igraph::degree(GBA)+1],
     vertex.size = 200 + (30 * igraph::degree(GBA)),
     vertex.label = NA)
par(oPar)                              # restore graphics state

# This is a very obviously different graph! Some biological networks have
# features that look like that - but in my experience the hub nodes are usually
# not that distinct. But then again, that really depends on the parameter
# "power". Play with the "power" parameter and get a sense for what difference
# this makes. Also: note that the graph has only a single component - no
# singletons.

# What's the degree distribution of this graph?
(dg <- igraph::degree(GBA))
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#DCF5B5",
     xlim = c(0,max(dg)+1), xaxt = "n",
     main = "Node degrees 200 nodes PA graph",
     xlab = "Degree", ylab = "Number")
axis(side = 1, at = seq(0, max(dg)+1, by=5))

# Most nodes have a degree of 1, but one node has a degree of 19.

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#DCF5B5",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a preferential-attachment network")

# Sort-of linear, but many of the higher ranked nodes have a frequency of only
# one. That behaviour smooths out in larger graphs:
#
X <- igraph::sample_pa(1e5, power = 0.8, directed = FALSE)  # 100,000 nodes
freqRank <- table(igraph::degree(X))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     xlab = "log(Rank)", ylab = "log(frequency)",
     pch = 21, bg = "#A5F5CC",
     main = "100,000 nodes in a random, scale-free network")
rm(X)


# ==   2.3  Random geometric graph  ============================================

# Finally, let's simulate a random geometric graph and look at the degree
# distribution. Remember: these graphs have a high probability to have edges
# between nodes that are "close" together - an entirely biological notion.

# We'll randomly place our nodes in a box. Then we'll define the
# probability for two nodes to have an edge to be a function of their Euclidian
# distance in the box.

# Here is a function that makes an adjacency matrix for such graphs. iGraph has
# a similar function, sample_grg(), which connects nodes that are closer than a
# cutoff, the function I give you below is a bit more interesting since it
# creates edges according to a probability that is determined by a generalized
# logistic function of the distance. This sigmoidal function gives a smooth
# cutoff and creates more "natural" graphs. Otherwise, the function is very
# similar to the random graph function, except that we output the "coordinates"
# of the nodes together with the adjacency matrix which we then use for the
# layout. list() FTW.
#

makeRandomGeometricAM <- function(nam, B = 25, Q = 0.001, t = 0.6) {
  # Make an adjacency matrix for an undirected random geometric graph from
  #    edges connected with probabilities according to a generalized logistic
  #    function.
  # Parameters:
  #    nam: a character vector of unique names
  #    B, Q, t:   probability that a random pair (i, j) of nodes gets an
  #                 edge determined by a generalized logistic function
  #                 p <- 1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / 0.9)))
  #
  # Value: a list with the following components:
  #        AM$mat : an adjacency matrix
  #        AM$nam : labels for the nodes
  #        AM$x   : x-coordinates for the nodes
  #        AM$y   : y-coordinates for the nodes
  #
  nu <- 1  # probably not useful to change
  AM <- list()
  AM$nam <- nam
  N <- length(AM$nam)
  AM$mat <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
  rownames(AM$mat) <- AM$nam
  colnames(AM$mat) <- AM$nam
  AM$x <- runif(N) # Randomly place nodes into the unit square
  AM$y <- runif(N)
  for (iRow in 1:(N-1)) { # Same principles as in makeRandomGraph()
    for (iCol in (iRow+1):N) {
      # geometric distance ...
      d <- sqrt((AM$x[iRow] - AM$x[iCol])^2 +
                (AM$y[iRow] - AM$y[iCol])^2)  # Pythagoras
      # distance dependent probability
      p <- 1 - 1/((1 + (Q * (exp(-B * (d-t)))))^(1 / nu))
      if (runif(1) < p) {
        AM$mat[iRow, iCol] <- 1
      }
    }
  }
  return(AM)
}

# Getting the parameters of a generalized logistic right takes a bit of
# experimenting. If you are interested, you can try a few variations. Or you can
# look up the function at
# https://en.wikipedia.org/wiki/Generalised_logistic_function

# This function computes generalized logistics ...
# genLog <- function(x, B = 25, Q = 0.001, t = 0.5) {
#     # generalized logistic (sigmoid)
#     nu <- 1
#     return(1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / nu)))
# }
#
# ... and this code plots p-values over the distances we could encouter between
# our nodes: from 0 to sqrt(2) i.e. the diagonal of the unit sqaure in which we
# will place our nodes.
# x <- seq(0, sqrt(2), length.out = 50)
# plot(x, genLog(x), type="l", col="#AA0000", ylim = c(0, 1),
#      xlab = "d", ylab = "p(edge)")

# 200 node random geomteric graph
N <- 200
set.seed(112358)                       # set RNG seed for repeatable randomness
rGAM <- makeRandomGeometricAM(as.character(1:N), t = 0.4)
set.seed(NULL)                         # reset the RNG


myGRG <- igraph::graph_from_adjacency_matrix(rGAM$mat, mode = "undirected")

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myGRG,
     layout = cbind(rGAM$x, rGAM$y), # use our node coordinates for layout,
     rescale = FALSE,
     xlim = c(min(rGAM$x) * 0.9, max(rGAM$x) * 1.1),
     ylim = c(min(rGAM$y) * 0.9, max(rGAM$y) * 1.1),
     vertex.color=heat.colors(max(igraph::degree(myGRG)+1))[igraph::degree(myGRG)+1],
     vertex.size = 0.1 + (0.2 * igraph::degree(myGRG)),
     vertex.label = NA)
par(oPar)

# degree distribution:
(dg <- igraph::degree(myGRG))
brk <- seq(min(dg) - 0.5, max(dg) + 0.5, by = 1)
hist(dg, breaks = brk, col = "#FCC6D2",
     xlim = c(0, 25), xaxt = "n",
     main = "Node degrees: 200 nodes RG graph",
     xlab = "Degree", ylab = "Number")
axis(side = 1, at = c(0, min(dg):max(dg)))

# You'll find that this is kind of in-between the random, and the scale-free
# graph. We do have hubs, but they are not as extreme as in the scale-free case;
# and we have have no singletons, in contrast to the random graph.

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#FCC6D2",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a random geometric network")



# =    3  A CLOSER LOOK AT THE igraph PACKAGE  =================================


# ==   3.1  Basics  ============================================================

# The basic object of the igraph package is a graph object. Let's explore the
# first graph some more, the one we built with our random gene names:
summary(myG)

# This output means: this is an IGRAPH graph, with U = UN-directed edges
#  and N = named nodes, that has 20 nodes and 20 edges. For details, see
?igraph::print.igraph

mode(myG)
class(myG)

# This means an igraph graph object is a special list object; it is opaque in
# the sense that a user is never expected to modify its components directly, but
# through a variety of helper functions which the package provides. There are
# many ways to construct graphs - from adjacency matrices, as we have just done,
# from edge lists, or by producing random graphs according to a variety of
# recipes, called _games_ in this package.

# Two basic functions retrieve nodes "Vertices", and "Edges":
igraph::V(myG)
igraph::E(myG)

# additional properties can be retrieved from the Vertices ...
igraph::V(myG)$name


# As with many R objects, loading the package provides special functions that
# can be accessed via the same name as the basic R functions, for example:

print(myG)
plot(myG)  # this is the result of default plot parameters

# ... where plot() allows the usual flexibility of fine-tuning the plot. We
# first layout the node coordinates with the Fruchtermann-Reingold algorithm - a
# force-directed layout that applies an ettractive potential along edges (which
# pulls nodes together) and a repulsive potential to nodes (so they don't
# overlap). Note the use of the degree() function to color and scale nodes and
# labels by degree and the use of the V() function to retrieve the vertex names.
# See ?plot.igraph for details."

# Plot with some customizing parameters
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myG,
     layout = igraph::layout_with_fr(myG),
     vertex.color=heat.colors(max(igraph::degree(myG)+1))[igraph::degree(myG)+1],
     vertex.size = 9 + (2 * igraph::degree(myG)),
     vertex.label.cex = 0.5 + (0.05 * igraph::degree(myG)),
     edge.width = 2,
     vertex.label = igraph::V(myG)$name,
     vertex.label.family = "sans",
     vertex.label.cex = 0.9)
par(oPar)

# ... or with a different layout:
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myG,
     layout = igraph::layout_in_circle(myG),
     vertex.color=heat.colors(max(igraph::degree(myG)+1))[igraph::degree(myG)+1],
     vertex.size = 9 + (2 * igraph::degree(myG)),
     vertex.label.cex = 0.5 + (0.05 * igraph::degree(myG)),
     edge.width = 2,
     vertex.label = igraph::V(myG)$name,
     vertex.label.family = "sans",
     vertex.label.cex = 0.9)
par(oPar)

# igraph has a large number of graph-layout functions: see
# ?layout_  and try them all.


# ==   3.2  Components  ========================================================

# The igraph function components() tells us whether there are components of the
# graph in which there is no path to other components.
igraph::components(myG)

# In the _membership_ vector, nodes are annotated with the index of the
# component they are part of. Sui7 is the only node of component 2, Cyj1 is in
# the third component etc. This is perhaps more clear if we sort by component
# index
sort(igraph::components(myG)$membership, decreasing = TRUE)

# Retrieving e.g. the members of the first component from the list can be done by subsetting:

(sel <- igraph::components(myG)$membership == 1)  # boolean vector ..
(c1 <- igraph::components(myG)$membership[sel])
names(c1)


# =    4  RANDOM GRAPHS AND GRAPH METRICS  =====================================


# Let's explore some of the more interesting, topological graph measures. We
# start by building a somewhat bigger graph. We aren't quite sure whether
# biological graphs are small-world, or random-geometric, or
# preferential-attachment ... but igraph has ways to simulate the basic ones
# (and we could easily simulate our own). Look at the following help pages:

?igraph::sample_gnm             # see also sample_gnp for the Erdös-Rényi models
?igraph::sample_smallworld      # for the Watts & Strogatz model
?igraph::sample_pa              # for the Barabasi-Albert model

# But note that there are many more sample_ functions. Check out the docs!

# Let's look at betweenness measures for our first graph. Here: the nodes again
# colored by degree. Degree centrality states: nodes of higher degree are
# considered to be more central. And that's also the way the force-directed
# layout drawas them, obviously.

set.seed(112358)                       # set RNG seed for repeatable randomness
myGxy <- igraph::layout_with_fr(myG)   # calculate layout coordinates
set.seed(NULL)                         # reset the RNG

oPar <- par(mar = rep(0, 4))           # turn margins off, save graphics state
plot(myG,
     layout = myGxy,
     rescale = FALSE,
     xlim = c(min(myGxy[,1]) * 0.99, max(myGxy[,1]) * 1.01),
     ylim = c(min(myGxy[,2]) * 0.99, max(myGxy[,2]) * 1.01),
     vertex.color=heat.colors(max(igraph::degree(myG)+1))[igraph::degree(myG)+1],
     vertex.size = 20 + (10 * igraph::degree(myG)),
     vertex.label = igraph::V(myG)$name,
     vertex.label.family = "sans",
     vertex.label.cex = 0.8)
par(oPar)                              # restore graphics state

# ==   4.1  Diameter  ==========================================================

igraph::diameter(myG)  # The diameter of a graph is its maximum length
                       # shortest path.

# let's plot this path: here are the nodes ...
igraph::get_diameter(myG)

# ... and we can get the x, y coordinates from iGxy by subsetting with the node
# names. The we draw the diameter-path with a transparent, thick pink line:
lines(myGxy[igraph::get_diameter(myG),], lwd=10, col="#ff63a788")

# == Centralization scores

?igraph::centralize
# replot our graph, and color by log_betweenness:

bC <- igraph::centr_betw(myG)  # calculate betweenness centrality
nodeBetw <- bC$res
nodeBetw <- round(log(nodeBetw +1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myG,
     layout = myGxy,
     rescale = FALSE,
     xlim = c(min(myGxy[,1]) * 0.99, max(myGxy[,1]) * 1.01),
     ylim = c(min(myGxy[,2]) * 0.99, max(myGxy[,2]) * 1.01),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 20 + (10 * igraph::degree(myG)),
     vertex.label = igraph::V(myG)$name,
     vertex.label.family = "sans",
     vertex.label.cex = 0.7)
par(oPar)

# Note that the betweenness - the number of shortest paths that pass through a
# node, is in general higher for high-degree nodes - but not always: Eqr2 has
# higher betweenness than Itv7: this measure really depends on the detailed
# local topology of the graph."

# Can you use centr_eigen() and centr_degree() to calculate the respective
# values? That's something I would expect you to be able to do.
#
# Lets plot betweenness centrality for our random geometric graph:

bCmyGRG <- igraph::centr_betw(myGRG)  # calculate betweenness centrality

nodeBetw <- bCmyGRG$res
nodeBetw <- round((log(nodeBetw +1))^2.5) + 1

# colours and size proportional to betweenness
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myGRG,
     layout = cbind(rGAM$x, rGAM$y), # use our node coordinates for layout,
     rescale = FALSE,
     xlim = c(min(rGAM$x) * 0.9, max(rGAM$x) * 1.1),
     ylim = c(min(rGAM$y) * 0.9, max(rGAM$y) * 1.1),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 0.1 + (0.03 * nodeBetw),
     vertex.label = NA)
par(oPar)

igraph::diameter(myGRG)
lines(rGAM$x[igraph::get_diameter(myGRG)],
      rGAM$y[igraph::get_diameter(myGRG)],
      lwd = 10,
      col = "#ff335533")



# =    5  GRAPH CLUSTERING  ====================================================


# Clustering finds "communities" in graphs - and depending what the edges
# represent, these could be complexes, pathways, biological systems or similar.
# There are many graph-clustering algorithms. One approach with many attractive
# properties is the Map Equation, developed by Martin Rosvall. See:
# http://www.ncbi.nlm.nih.gov/pubmed/18216267 and htttp://www.mapequation.org


myGRGclusters <- igraph::cluster_infomap(myGRG)
igraph::modularity(myGRGclusters)  # ... measures how separated the different
                                   # membership types are from each other
igraph::membership(myGRGclusters)         # which nodes are in what cluster?
table(igraph::membership(myGRGclusters))  # how large are the clusters?

# The largest cluster has 48 members, the second largest has 25, etc.


# Lets plot our graph again, coloring the nodes of the first five communities by
# their cluster membership:

# first, make a vector with as many grey colors as we have communities ...
commColors <- rep("#f1eef6", max(igraph::membership(myGRGclusters)))
# ... then overwrite the first five with "real colors" - something like rust,
# lilac, pink, and mauve or so.
commColors[1:5] <- c("#980043", "#dd1c77", "#df65b0", "#c994c7", "#d4b9da")


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(myGRG,
     layout = cbind(rGAM$x, rGAM$y),
     rescale = FALSE,
     xlim = c(min(rGAM$x) * 0.9, max(rGAM$x) * 1.1),
     ylim = c(min(rGAM$y) * 0.9, max(rGAM$y) * 1.1),
     vertex.color=commColors[igraph::membership(myGRGclusters)],
     vertex.size = 0.1 + (0.1 * igraph::degree(myGRG)),
     vertex.label = NA)
par(oPar)




# [END]
