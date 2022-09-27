# geneLengthAcrossOrganisms.R
#
# Code from class session 02, 2022-09-20


# Task: plat the gene length distribution for a model organism.
# Pseudocode:
# ===========
#
# Access the genome
# Retrieve annotations
# Retrieve start and end coordinates of the gene
# Compute the length and store data
# Plot data as histogram
#
# What are possible sources for gene annotations?
#
#   - Bioconductor goseq:: package has  getlength() function ... but
#       this requires a vector of genes as input. Not useful for our purposes.
#
#   - Genbank ? Should have it. We couldn't find this information easily.
#        N.b.  On the NCBI homepage (https://www.ncbi.nlm.nih.gov/)
#              you can find a link to "Genomes & Maps". A link to the Genome
#              takes you to https://www.ncbi.nlm.nih.gov/genome . There, a
#              search by organism finds an overview page for the yeast genome.
#                  https://www.ncbi.nlm.nih.gov/genome/15
#              This contains a lot of information, which we should have looked
#              at systematically: there is indeed a line that says: "Download
#              genome annotation in GFF, GenBank or tabular format" ... and each
#              of these formats contain what we were looking for.
#
#              Bottom line: Genbank is one of two authoritative databases for
#              genome information, though you need to search carefully.
#
#   - UCSC genome browser?
#     UCSC is another great genome annotation site. At the home page
#        https://genome.ucsc.edu/
#     choose the Genomes > Other link and search for
#     "Saccharomyces cerevisiae". Careful reading of the page finds (towards
#     the bottom): "Downloads of the yeast data and annotations may be obtained
#     from the UCSC Genome Browser FTP server or Downloads page. "
#
#     The Downloads page (http://hgdownload.soe.ucsc.edu/downloads.html#yeast)
#     has several annotation files. But it is easy to find which ones are
#     useful to our quest. Just Google for "2bit genome annotation format",
#     "GTF genome annotation format" etc. and it should quickly become obvious.
#
#     Ultimately we downloaded a GTF formatted annotation file. It had the
#     extension  .gtf.gz  ... and ".gz" indicates one of the widely used
#     file-compression formats.
#     On my Mac, I just double click .gz files to unpack them.
#     On Windows?
#     On linux I do something like    $ gunzip infile.gz
#
#
#   - EBI? The European Bioinformatics Institute is the second authoritative
#     source of molecular data.  https://www.ebi.ac.uk/
#
#     Typing "genome" into their database and tool finding search field finds
#     the very excellent ensembl! resource. http://useast.ensembl.org/index.html
#     (This is a North-American mirror.) This has a link to BioMart
#       http://useast.ensembl.org/biomart/martview/
#     I versatile web form interface allows you to produce custom queries and
#     format custom output files, and within a few minutes you can produce
#     the following data and download it (slightly edited):
#
#     Gene ID	   Start     End      Chr.   Name
#     YBR024W    289445    290350	  II     SCO2
#     YDL245C    11657     13360    IV     HXT15
#     YBR021W    281443    283344	  II     FUR4
#     YGR014W    516943    520863	  VII    MSB2
#     YKL119C    218570    219217	  XI     VPH2
#     YPR031W    631515    633761	  XVI    NTO1
#     YJL076W    295245    298814	  X      NET1
#     YPL174C    220167    222773	  XVI    NIP100
#     [... 5805 rows]
#
#     That would have been my favourite data to work with.

#
# We downloaded GCF_000146045.2_R64.augustus.gtf from UCSC, I
# placed it into ./data/GCF_000146045.2_R64.augustus.gtf and used a standard
# R invocation to read it.

myGTF <- read.delim("./data/GCF_000146045.2_R64.augustus.gtf",
                    header = FALSE)

# The variable myGTF is a "dataframe" - one of the most frequently used
# R datastructures, essentially a table of data.
#
head(myGTF)

# Column V3 seems to have the name of the feature. How many different features
# are there?
table(myGTF$V3)

#   CDS    start_codon     stop_codon
#   5353          4964           4963

# We can filter myGTF to contain only those rows which have the string "CDS"
# in column V3.
myCDS <- myGTF[myGTF$V3 == "CDS", ]

head(myCDS)   # confirm...

# Now, compute the gene lengths ...
# test: myCDS$V5[1:20] - myCDS$V4[1:20]
myLengths <- myCDS$V5 - myCDS$V4

# ... analyze the results
summary(myLengths)

# ... and plot a histogram.
hist(log10(myLengths))

# This is the gene-length distribution we were looking for.
# ---------------------------------------------------------
#
#
# PS: Just for fun we also looked at the number of genes in each chromosome.
barplot(sort(table(myCDS$V1)))


# [END]
