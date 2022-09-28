# geneLengthFromBiomart.R
#
# Code from class session 03, 2022-09-27

# ==============================================================================

# Task: Use the ensembl biomart web interface to download gene data. Load it
# in R, check integrity, plot

# - Access http://ebi.ac.uk
#
# - Many services ... follow "Find data resources", search for "genome".
#
# - The top hit - ensembl - takes us to the ensembl genome database, with its
#   API and data resources. The link redirects you to a mirror site (Nb. we
#   usually get taken to "useast", but there may have been a service disruption
#   since we were taken to "uswest" and the response times were slow.)
#
# - We followed a link to Tools > Biomart (The text that looked promising was:
#   "Export custom datasets...") http://uswest.ensembl.org/biomart/martview/
#
# - This is the top entry point for a forms-based Web interface to the ensmbl
#   genomes. From here, we need to make the following choices to download data:
#   - Dataset (which genome)
#   - Filters (which genes from that genome)
#   - Attributes (which features of each gene)
#
#   We selected:
#   - Dataset: Ensembl genes -> Saccaromyces cerevisiae
#   - Filters: Gene -> Limit to genes (external references)... -> With SGD
#              gene name IDs  (SGD is Saccharomyces Genome Database
#                              https://www.yeastgenome.org)
#   - Attributes: Gene -> Gene stable ID
#                         Transcript stable ID
#                         Gene description
#                         Gene start (bp)
#                         Gene end (bp)
#                         Strand
#                         Transcript start (bp)
#                         Transcript end (bp)
#                         Gene name
#                       * Transcript count
#                         Gene type
#
#                         (*) We were not sure what this is, and could not find
#                             documentation how this is defined. Nb. many data
#                             items are defined at
#                       http://ensembl.org/info/docs/api/core/core_schema.html
#                             ... but not this one.

# We first examined the resulting data set on its Web page, then downloaded the
# data in csv (comma separated values) format. It is stored in
#  ./data/S.cerevisiae.GenomeGenesBiomart.csv
#
#
# == 1: Load the data ==========================================================

# Pseudocode:
# ===========
#
# load data as dataframe
SCgeneData <- read.csv("./data/S.cerevisiae.GenomeGenesBiomart.csv")
View(SCgeneData)


# == 2: Integrity checks =======================================================

# integrity checks:
#    - is the start coordinate always smaller than the end coordinate?
#    subtract start from end in each row
x <- SCgeneData$Transcript.start..bp.         # (two alternative but equivalent
y <- SCgeneData[ , "Transcript.end..bp."]     #  idioms to access data from a
                                              # data frame)
head(x)
head(y)
all(x < y)

#    - are there unobserved values?
for (i in 1:ncol(SCgeneData)) {
  print(any(is.na(SCgeneData[ , i])))
}

#    - what transcript types do we see? are there pseudogenes?
unique(SCgeneData$Gene.type)
table(SCgeneData$Gene.type)

#    - are they all yeast genes?
x <- ! grepl("^Y", SCgeneData$Gene.stable.ID)
sum(x)

SCgeneData$Gene.stable.ID[x]

#    - are there very short genes? Suspiciously long genes?
#      - compute and store the lengths
x <- SCgeneData$Transcript.start..bp.
y <- SCgeneData$Transcript.end..bp.
SCgeneData$length <- y - x + 1   # add  one to difference to get length. Note
                                 # that R automatically creates a new column
                                 # when we assign values to it.

plot(sort(log(SCgeneData$length)))

hist(log10(SCgeneData$length), breaks=20)

summary(SCgeneData$length)

SCgeneData[SCgeneData$length == min(SCgeneData$length) , ]
SCgeneData[SCgeneData$length == max(SCgeneData$length) , ]



# == 3: More .. ? ==============================================================

# - function to do this automatically given a species name as input
# - function to ask about functional annotations as a function of gene length
# - predict gene category from length?
#



# [END]
