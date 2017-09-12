# FND-MAT-Graphs_and_networks.R
#
# Purpose:  A Bioinformatics Course:
#              R code accompanying the FND-MAT-Graphs_and_networks unit.
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

# This tutorial covers basic concepts of graph theory and analysis in R. You
# should have typed init() to configure some utilities in the background.


# ==============================================================================
#        PART ONE: REVIEW
# ==============================================================================

# I assume you'll have read the Pavlopoulos review of graph theory concepts.
# Let's explore some of the ideas by starting with a small random graph."


# To begin let's write a little function that will create random "gene" names;
# there's no particular purpose to this other than to make our graphs look a
# little more like what we would find in a publication ...
makeRandomGenenames <- function(N) {
  nam <- character()
  while (length(nam) < N) {
    a <- paste(c(sample(LETTERS, 1), sample(letters, 2)),
               sep="", collapse="") # three letters
    n <- sample(1:9, 1)             # one number
    nam[length(nam) + 1] <- paste(a, n, sep="") # store in vector
    nam <- unique(nam)   # delete if this was a duplicate
  }
  return(nam)
}

N <- 20

set.seed(112358)
Nnames <- makeRandomGenenames(N)

Nnames

# One way to represent graphs in a computer is as an "adjacency matrix". In this
# matrix, each row and each column represents a node, and the cell at the
# intersection of a row and column contains a value/TRUE if there is an edge,
# 0/FALSE otherwise. It's easy to see that an undirected graph has a symmetric
# adjacency matrix (i, j) == (j, i); and we can put values other than {1, 0}
# into a cell if we want to represent a weighted edge.

# At first, lets create a random graph: let's say a pair of nodes has
# probability p <- 0.1 to have an edge, and our graph is symmetric and has no
# self-edges. We use our Nnames as node labels, but I've written the function so
# that we could also just ask for any number of un-named nodes, we'll use that later.

makeRandomGraph <- function(nam, p = 0.1) {
  # nam: either a character vector of unique names, or a single
  #        number that will be converted into a vector of integers.
  # p:   probability that a random pair of nodes will have an edge.
  #
  # Value: an adjacency matrix
  #
  if (is.numeric(nam) && length(nam) == 1) { # if nam is  a single number ...
    nam <- as.character(1:nam)
  }
  N <- length(nam)
  G <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
  rownames(G) <- nam
  colnames(G) <- nam
  for (iRow in 1:(N-1)) { # Note how we make sure iRow != iCol
    for (iCol in (iRow+1):N) {
      if (runif(1) < p) {  # runif() creates uniform random numbers
        # between 0 and 1
        G[iRow, iCol] <- 1   # row, col !
        G[iCol, iRow] <- 1   # col, row !
      }
    }
  }
  return(G)
}

set.seed(112358)
G <- makeRandomGraph(Nnames, p = 0.09)
G


# Listing the matrix is not very informative - we should plot this graph. We'll
# go into more details of the igraph package a bit later, for now we just use it
# to plot:

if (!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}

iG <- graph_from_adjacency_matrix(G)
iGxy <- layout_with_graphopt(iG, charge=0.001)   # calculate layout coordinates


# The igraph package adds its own function to the collection of plot()
# functions; R makes the selection which plot function to use based on the class
# of the object that we request to plot. This plot function has parameters
#  layout - the x,y coordinates of the nodes;
#  vertex.color - which I define to color by node-degree
#  vertex size - which I define to increase with node-degree
#  vertex.label - which I set to use our Nnames vector

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 800 + (150 * degree(iG)),
     vertex.label = as.character(degree(iG)/2),
     #     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)  # reset plot window


# The simplest descriptor of a graph are the number of nodes, edges, and the
# degree-distribution. In our example, the number of nodes was given: N; the
# number of edges can easily be calculated from the adjacency matrix. In our
# matrix, we have entered 1 for every edge. Thus we simply sum over the matrix:
sum(G)

# Is that correct? Is that what you see in the plot?

# Yes and no: we entered every edge twice: once for a node [i,j], and again for
# the node [j, i]. Whether that is correct depends on what exactly we
# want to do with the matrix. If these were directed edges, we would need to
# keep track of them separately. Since we didn't intend them to be directed,
# we'll could divide the number of edges by 2. Why didn't we simply use an
# upper-triangular matrix? Because then we need to keep track of the ordering of
# edges if we want to know whether a particular edge exists or not. For example
# we could sort the nodes alphabetically, and make sure we always query a pair
# in alphabetical order. Then a triangular matrix would be efficient.

# What about the degree distribution? We can get that simply by summing over the
# rows (or the columns):"

rowSums(G)  # check this against the plot!

# Let's  plot the degree distribution in a histogram:
rs <- rowSums(G)
brk <- seq(min(rs)-0.5, max(rs)+0.5, by=1)  # define breaks for the histogram
hist(rs, breaks=brk, col="#A5CCF5",
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

# ==============================================================================
#        PART TWO: DEGREE DISTRIBUTIONS
# ==============================================================================

# Let's simulate a few graphs that are a bit bigger to get a better sense of
# their degree distributions:
#

# === random graph


set.seed(31415927)
G200 <- makeRandomGraph(200, p = 0.015)
iG200 <- graph_from_adjacency_matrix(G200)
iGxy <- layout_with_graphopt(iG200, charge=0.0001) # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG200,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG200)+1))[degree(iG200)+1],
     vertex.size = 200 + (30 * degree(iG200)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# This graph has thirteen singletons and one large, connected component. Many
# biological graphs look approximately like this.

# Calculate degree distributions
dg <- degree(iG200)/2   # here, we use the iGraph function degree()
# not rowsums() from base R.
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5CCF5",
     xlim = c(-1,11), xaxt = "n",
     main = "Node degrees", xlab = "Degree", ylab = "Number")  # plot histogram
axis(side = 1, at = 0:10)



# Note the characteristic peak of this distribution: this is not "scale-free". Here is a log-log plot of frequency vs. degree-rank:

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#A5CCF5",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a random network")

# === scale-free graph (Barabasi-Albert)

# What does one of those intriguing "scale-free" distributions look like? The
# iGraph package has a function to make random graphs according to the
# Barabasi-Albert model of scale-free graphs. It is: sample_pa(), where pa
# stands for "preferential attachment", one type of process that will yield
# scale-free distributions.


set.seed(31415927)
GBA <- sample_pa(200, power = 0.8)

iGxy <- layout_with_graphopt(GBA, charge=0.0001) # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(GBA,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(GBA)+1))[degree(GBA)+1],
     vertex.size = 200 + (30 * degree(GBA)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# This is a very obviously different graph! Some biological networks have
# features that look like that - but in my experience the hub nodes are usually
# not that distinct. But then again, that really depends on the parameter
# "power". Feel encouraged to change "power" and get a sense for what difference
# this makes. Also: note that the graph has only a single component.

# What's the degree distribution of this graph?
(dg <- degree(GBA))
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5D5CC",
     xlim = c(0,30), xaxt = "n",
     main = "Node degrees 200 nodes PA graph",
     xlab = "Degree", ylab = "Number")
axis(side = 1, at = seq(0, 30, by=5))

# Most nodes have a degree of 1, but one node has a degree of 28.

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#A5F5CC",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a preferential-attachment network")

# Sort-of linear, but many of the higher ranked nodes have a frequency of only
# one. That behaviour smooths out in larger graphs:
#
X <- sample_pa(100000, power = 0.8)  # 100,000 nodes
freqRank <- table(degree(X))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     xlab = "log(Rank)", ylab = "log(frequency)",
     pch = 21, bg = "#A5F5CC",
     main = "100,000 nodes in a random, scale-free network")
rm(X)

# === Random geometric graph

# Finally, let's simulate a random geometric graph and look at the degree
# distribution. Remember: these graphs have a high probability to have edges
# between nodes that are "close" together - an entriely biological notion.

# We'll randomly place our nodes in a box. Then we'll define the
# probability for two nodes to have an edge to be a function of their distance.

# Here is a function that makes such graphs. iGraph has sample_grg(), which
# connects nodes that are closer than a cutoff, the function I give you below is
# a bit more interesting since it creates edges according to a probability that
# is determined by a generalized logistic function of the distance. This
# sigmoidal function gives a smooth cutoff and creates more "natural" graphs.
# Otherwise, the function is very similar to the random graph function, except
# that we output the "coordinates" of the nodes together with the adjacency
# matrix. Lists FTW.
#
makeRandomGeometricGraph <- function(nam, B = 25, Q = 0.001, t = 0.6) {
  # nam: either a character vector of unique names, or a single
  #        number that will be converted into a vector of integers.
  # B, Q, t:   probability that a random pair (i, j) of nodes gets an
  #              edge determined by a generalized logistic function
  #              p <- 1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / 0.9)))
  #
  # Value: a list with the following components:
  #        G$mat : an adjacency matrix
  #        G$nam : labels for the nodes
  #        G$x   : x-coordinates for the nodes
  #        G$y   : y-coordinates for the nodes
  #
  nu <- 1  # probably not useful to change
  G <- list()

  if (is.numeric(nam) && length(nam) == 1) {
    nam <- as.character(1:nam)
  }
  G$nam <- nam
  N <- length(G$nam)
  G$mat <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
  rownames(G$mat) <- G$nam
  colnames(G$mat) <- G$nam
  G$x <- runif(N)
  G$y <- runif(N)
  for (iRow in 1:(N-1)) { # Same principles as in makeRandomGraph()
    for (iCol in (iRow+1):N) {
      # geometric distance ...
      d <- sqrt((G$x[iRow] - G$x[iCol])^2 +
                  (G$y[iRow] - G$y[iCol])^2)  # Pythagoras
      # distance dependent probability
      p <- 1 - 1/((1 + (Q * (exp(-B * (d-t)))))^(1 / nu))
      if (runif(1) < p) {
        G$mat[iRow, iCol] <- 1
        G$mat[iCol, iRow] <- 1
      }
    }
  }
  return(G)
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
set.seed(112358)
GRG <- makeRandomGeometricGraph(200, t=0.4)


iGRG <- graph_from_adjacency_matrix(GRG$mat)
iGRGxy <- cbind(GRG$x, GRG$y) # use our node coordinates for layout

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])) * 1.1,
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iGRG)+1))[degree(iGRG)+1],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# degree distribution:
(dg <- degree(iGRG)/2)
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#FCD6E2",
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
     pch = 21, bg = "#FCD6E2",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a random geometric network")



# ====================================================================
#        PART THREE: A CLOSER LOOK AT THE igraph PACKAGE
# ====================================================================


# == BASICS ==========================================================

# The basic object of the igraph package is a graph object. Let's explore the
# first graph some more, the one we built with our random gene names:
summary(iG)

# This output means: this is an IGRAPH graph, with D = directed edges and N =
# named nodes, that has 20 nodes and 40 edges. For details, see
?print.igraph

mode(iG)
class(iG)

# This means an igraph graph object is a special list object; it is opaque in
# the sense that a user is never expected to modify its components directly, but
# through a variety of helper functions which the package provides. There are
# many ways to construct graphs - from adjacency matrices, as we have just done,
# from edge lists, or by producing random graphs according to a variety of
# recipes, called _games_ in this package.

# Two basic functions retrieve nodes "Vertices", and "Edges":
V(iG)
E(iG)

# As with many R objects, loading the package provides special functions that
# can be accessed via the same name as the basic R functions, for example:

print(iG)
plot(iG)

# ... where plot() allows the usual flexibility of fine-tuning the plot. We
# first layout the node coordinates with the Fruchtermann-Reingold algorithm - a
# force-directed layout that applies an ettractive potential along edges (which
# pulls nodes together) and a repulsive potential to nodes (so they don't
# overlap). Note the use of the degree() function to color and scale nodes and
# labels by degree and the use of the V() function to retrieve the vertex names.
# See ?plot.igraph for details."

iGxy <- layout_with_fr(iG)   # calculate layout coordinates

# Plot with some customizing parameters
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 9 + (2 * degree(iG)),
     vertex.label.cex = 0.5 + (0.05 * degree(iG)),
     edge.arrow.size = 0,
     edge.width = 2,
     vertex.label = toupper(V(iG)$name))
par(oPar)


# == Components

# The igraph function components() tells us whether there are components of the
# graph in which there is no path to other components.
components(iG)

# In the _membership_ vector, nodes are annotatd with the index of the component
# they are part of. Sui7 is the only node of component 2, Cyj1 is in the third
# component etc. This is perhaps more clear if we sort by component index
sort(components(iG)$membership)

# Retrieving e.g. the members of the first component from the list can be done by subsetting:

components(iG)$membership == 1  # logical ..
components(iG)$membership[components(iG)$membership == 1]
names(components(iG)$membership)[components(iG)$membership == 1]



# == RANDOM GRAPHS AND GRAPH METRICS =================================


# Let's explore some of the more interesting, topological graph measures. We
# start by building a somewhat bigger graph. We aren't quite sure whether
# biological graphs are small-world, or random-geometric, or
# preferential-attachment ... but igraph has ways to simulate the basic ones
# (and we could easily simulate our own). Look at the following help pages:

?sample_gnm                      # see also sample_gnp for the Erdös-Rényi models
?sample_smallworld               # for the Watts & Strogatz model
?sample_pa                       # for the Barabasi-Albert model

# But note that there are many more sample_ functions. Check out the docs!

# Let's look at betweenness measures for our first graph: here: the nodes again
# colored by degree. Degree centrality states: nodes of higher degree are
# considered to be more central. And that's also the way the force-directed
# layout drawas them, obviously.

set.seed(112358)
iGxy <- layout_with_fr(iG)   # calculate layout coordinates
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 20 + (10 * degree(iG)),
     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)

# == Diameter

diameter(iG)  # The diameter of a graph is its maximum length shortest path.

# let's plot this path: here are the nodes ...
get_diameter(iG)

# ... and we can get the x, y coordinates from iGxy by subsetting with the node
# names. The we draw the diameter-path with a transparent, thick pink line:
lines(iGxy[get_diameter(iG),], lwd=10, col="#ff63a788")

# == Centralization scores

?centralize
# replot our graph, and color by log_betweenness:

bC <- centr_betw(iG)  # calculate betweenness centrality
nodeBetw <- bC$res
nodeBetw <- round(log(nodeBetw +1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 20 + (10 * degree(iG)),
     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)

# Note that the betweenness - the number of shortest paths that pass through a
# node, is in general higher for high-degree nodes - but not always: Eqr2 has
# higher betweenness than Itv7: this measure really depends on the detailed
# local topology of the graph."

# Can you use centr_eigen() and centr_degree() to calculate the respective
# values? That's something I would expect you to be able to do.
#
# Lets plot betweenness centrality for our random geometric graph:

bCiGRG <- centr_betw(iGRG)  # calculate betweenness centrality

nodeBetw <- bCiGRG$res
nodeBetw <- round((log(nodeBetw +1))^2.5) + 1

# colours and size proportional to betweenness

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])),
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 0.1 + (0.03 * nodeBetw),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

diameter(iGRG)
lines(iGRGxy[get_diameter(iGRG),], lwd=10, col="#ff335533")



# == CLUSTERING ======================================================

# Clustering finds "communities" in graphs - and depending what the edges
# represent, these could be complexes, pathways, biological systems or similar.
# There are many graph-clustering algorithms. One approach with many attractive
# properties is the Map Equation, developed by Martin Rosvall. See:
# http://www.ncbi.nlm.nih.gov/pubmed/18216267 and htttp://www.mapequation.org


iGRGclusters <- cluster_infomap(iGRG)
modularity(iGRGclusters) # ... measures how separated the different membership
# types are from each other
membership(iGRGclusters) # which nodes are in what cluster?
table(membership(iGRGclusters))  # how large are the clusters?

# The largest cluster has 48 members, the second largest has 25, etc.


# Lets plot our graph again, coloring the nodes of the first five communities by
# their cluster membership:

# first, make a vector with as many grey colors as we have communities ...
commColors <- rep("#f1eef6", max(membership(iGRGclusters)))
# ... then overwrite the first five with "real colors" - something like rust,
# lilac, pink, and mauve or so.
commColors[1:5] <- c("#980043", "#dd1c77", "#df65b0", "#c994c7", "#d4b9da")


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])),
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])),
     vertex.color=commColors[membership(iGRGclusters)],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)


# = 1 Tasks




# [END]
