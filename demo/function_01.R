# function_01.R
#
#   Starting code in class session 07, 2022-10-25
#
#   Analysis motivated by Abou Chakra et al. (2021)
#    Control of tissue development and cell diversity by
#    cell cycle-dependent transcriptional filtering
#      https://elifesciences.org/articles/64951
#

# ==============================================================================

# Task:
#    I: Definitions of function
#   II: What data can we obtain
#  III: How do we find "function" in the data

#
# ==============================================================================
#    I: What is function?
# ==============================================================================
#

# In Figure 8, Abou-Chakra et al. suggest that short human genes have different
# functional annotations enriched than long genes, suggesting a functional role for transcriptional filtering. One could imagine that genes that are expressed at later time points in development (longer cell-cycle times) ought to be systematically different from genes that are expressed in earlier developmental stages, but does transcriptional filtering really contribute significantly to such regulation, and what functions exactly are involved in this?

# We need to clarify:
#   - What is "function"? How can we obtain a computable concept of function?
#   - What do we mean by enrichment?


# == 1.1 Function ...

# Every gene has a function. At least one. We know thatm because all genes
# undergo random mutations, and it is only purifying selection that retains the
# gene. We usually think of purifying selection to be driven by a gene's
# function, although things may be a lot more complicated in reality. But this
# is a valid starting assumption.

# However, that does not allow us to ask a question like: "do short genes have
# different functions from long genes?" The question makes no sense - a gene has
# the function it has and the length it has and that's it. All functions are
# different. If they were not, genes would be redundant and purifying selection
# would not work. In  order for the question to be meaningful, we need to
# establish groups of functions that are similar. About-Chakra et al. call these
# groupings "themes". And if you can define a functional theme, you indeed can
# ask: "do genes that are annotated to a theme have a higher likelihood of being
# short or long, relative to their overall occurrence in the genome?" Or: "Are
# there themes that are enriched among short genes?" (Which is not exactly the
# same question.)

# How are themes defined, and how do we establish which genes are annotated to a
# theme?

# == 1.1 Enrichment ...

#
# ==============================================================================
#    II: Data
# ==============================================================================
#

inFile <- "demo/elife-64951-supp5-v2.csv"
file.exists(inFile)
dat1 <- read.csv(inFile)

any(duplicated(dat1$Theme))  # False
x <- unlist(dat1[ , -1])
any(duplicated(x))  # False
 myTab <- sort(table(x), decreasing = TRUE)
head(myTab)
s <- names(myTab)[1]
x <- x[x != ""]
myTab <- sort(table(x), decreasing = TRUE)
head(myTab)

for (i in seq_len(nrow(dat1))) {
  idx <- grep("cardiocyte", dat1[i, ])
  if (length(idx) > 0) {
    print(sprintf("found in row: %d", i))
  }
}

dat1[c(9,16,29), ]

grep("cardiocyte", dat1[i, 9])

length(myTab[myTab > 1])

# =============
#
inFile <- "demo/elife-64951-supp6-v2.csv"
file.exists(inFile)
dat2 <- read.csv(inFile)

# ====
# Next: read this into a list
#       define the semantics of dat2
#       see how to crossreference this with actual genes
#       figure out the process of associating an actual gene with
#       its function theme(s)
#       look for enrichment
#       produce a wordcloud


# [END]
