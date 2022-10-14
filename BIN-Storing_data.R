# tocID <- "BIN-Storing_data.R"
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-Storing_data unit
#
# Version: 1.4.2
#
# Date:    2017-10  -  2022-10
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.4.2  fixed broken UniProt ID mapping link (thanks Tatsuya Corlett)
# V 1.4.1  removed a stray, forgotten formatting instruction
# V 1.4    2022 removed "submit for credit"
# V 1.3.2  2021 minimal maintenance
# V 1.3.1  add overlooked  jsonlite:: prefix to fromJson()
# V 1.3    Made file locations more consistent. All student-edited files
#          go into the myScripts directory
# V 1.2    2020 updates. Finally removed stringAsFactors  :-)
# V 1.1    Add instructions to retrieve UniProt ID from ID mapping service.
# V 1.0    First live version, complete rebuilt. Now using JSON data sources.
# V 0.1    First code copied from BCH441_A03_makeYFOlist.R
#
# TODO:
#  The sameSpecies() approach is a bit of a hack - can we solve the
#  species vs. strain issue in a more principled way?
#
# == HOW TO WORK WITH LEARNING UNIT FILES ======================================
#
# DO NOT SIMPLY  source()  THESE FILES!
#
# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
#  going on. That's not how it works ...
#
# ==============================================================================


#TOC> ==========================================================================
#TOC>
#TOC>   Section  Title                                                   Line
#TOC> -----------------------------------------------------------------------
#TOC>   1        A Relational Datamodel in R: review                       66
#TOC>   1.1        Building a sample database structure                   106
#TOC>   1.1.1          completing the database                            212
#TOC>   1.2        Querying the database                                  245
#TOC>   1.3        Task: Exercise                                         276
#TOC>   2        Implementing the protein datamodel                       295
#TOC>   2.1        JSON formatted source data                             323
#TOC>   2.2        "Sanitizing" sequence data                             364
#TOC>   2.3        Create a protein table for our data model              388
#TOC>   2.3.1          Initialize the database                            390
#TOC>   2.3.2          Add data                                           402
#TOC>   2.4        Complete the database                                  422
#TOC>   2.4.1          Examples of navigating the database                449
#TOC>   2.5        Updating the database                                  481
#TOC>   3        Add your own data                                        493
#TOC>   3.1        Find a protein                                         501
#TOC>   3.2        Put the information into JSON files                    533
#TOC>   3.3        Create an R script to create your own database         575
#TOC>   3.3.1          Check and validate                                 604
#TOC>   3.4        Record the results in your journal:                    648
#TOC>
#TOC> ==========================================================================


# =    1  A Relational Datamodel in R: review  =================================

# A disclaimer at first: we are not building an industry-strength database at
# all here - but we are employing principles of such a database to keep common
# types of lab-data well organized. Don't think of this as emulating or even
# replacing a "real" database, but think of this as improving the ad hoc
# approaches we normally employ to store data in the lab. That does not mean
# such ad hoc approaches are necessarily bad - the best solution always depends
# on your objectives, the details of your tasks, and the context in which you
# are working.

# The principle we follow in implementing a relational data model is to build a
# list of dataframes . This list is our "database":
#  - Each _entity_ of the datamodel is a dataframe. In an SQL database, these
#    would also be called "tables". In a spreadsheet this would be a "sheet".
#  - Each instance of an entity, i.e. one stored _item_, is a row of the data
#    frame. In an SQL database this would be a record. In a spreadsheet this is
#    a row.
#  - Each _attribute_ of an entity is is a column of the dataframe. In an SQL
#    database this is a column, in a spreadsheet too.
#  - This doesn't necessarily solve the question of how we will store and curate
#    our source data - we will defer that to later. At first we talk only about
#    data representation internal to our R session, where we need it for
#    processing and analysis.

# Lets review syntax for creating and accessing such a structure, a list of data
# frames. You'll have to be absolutely confident with this, or you'll get lost
# in all the later learning units. We'll start from a compact example, a tiny
# database of philosophers to keep things brief. That database will have three
# tables: person, works and book. Person stores biographical data, book stores
# books, and works is a join table associating persons with their work. You
# should already be familiar with "join tables" and why we need them. This is
# the structure:
#
# person: id, name, born, died, school
# book:   id, title, published
# works:  id, person$id, book$id

# Perhaps draw out this schema to make things more clear.

# ==   1.1  Building a sample database structure  ==============================

# Let's build this structure.

philDB <- list()  # This is an empty list

# This is a data frame that we initialize with two philosophers
x <- data.frame(id = c(1,2),
                name = c("Laozi", "Martin Heidegger"),
                born = c(NA, "1889"),
                died = c("531 BCE", "1976"),
                school = c("Daoism", "Phenomenology"))
str(x)

# Lets add the dataframe to the philDB list and call it "person" there.
philDB[["person"]] <- x
str(philDB)

# and let's remove x so we don't mix up things later.
rm(x)

# We can address elements with the usual subsetting operators. I will use
# the $ operator for tables and columns, the [] operator for elements in
# columns. For example ...

philDB$person$name[1]   # Laozi

# task: Write an expression that returns all "school" entries from the
#       person table.

# Let's now add another person. There are several ways to do this, the
# conceptually cleanest is to create a one-row dataframe with the data, and
# rbind() it to the existing dataframe. Doing this, we must take care that
# the data frame column names are identical. What happens if they are not?
# Let's find out:
(x <- data.frame(a=1:4, b=11:14))
(y <- data.frame(a=6, c=17))
rbind(x, y)

(y <- data.frame(a=6, b=17))
rbind(x, y)

# All clear? That's good - this behaviour provides us with a sanity check on the
# operation. Incidentally: rbind(x, y) did NOT change the table ...
x
# rather rbind() had the changed table as its return value and that's why it
# was printed. To actually change the table, you need to ASSIGN the return
# value  of rbind() ... like so:
x <- rbind(x, y)

# To continue ...
(x <- data.frame(id = 2,
                 name = "Zhuangzi",
                 born = "369 BCE",
                 died = "286 BCE",
                 school = "Daoism"))

# Add this to the "person" table in our database with rbind() ...
philDB$person <- rbind(philDB$person, x)

# ... and examine the result:
str(philDB)

# We made a serious error in our data! Did you spot it?
#
# If not, look at ...
philDB$person$id
# ... does that look oK?
#
# Absolutely not! "id" is the Primary Key in the table, and it has to be
# unique. How can we guarantee it to be unique? Certainly not when we
# enter it by hand. We need a function that generates a unique key. Here's
# a simple version, without any error-checking. It assumes that a column
# named "id" exists in the table, and that it holds the Primary Keys:
autoincrement <- function(table) {
  return(max(table$id) + 1)
}

#Try it:
autoincrement(philDB$person)

# Once that is clear, let's remove the Zhuangzi entry and recreate it correctly.
# Many ways to remove, here we use a logical expression to select matching
# record(s), apply the results to subset the data frame, and overwrite the
# existing table with the new one.

sel <- !(philDB$person$name == "Zhuangzi")   # select ...
philDB$person <- philDB$person[sel, ]        # ... and replace

str(philDB)

# Now let's add Zhuangzi with correct data. Note how we use the autoincrement
# function for the id

x <- data.frame(id = autoincrement(philDB$person),
                name = "Zhuangzi",
                born = "369 BCE",
                died = "286 BCE",
                school = "Daoism")
philDB$person <- rbind(philDB$person, x)
str(philDB)

# So far so good. Be honest with yourself. If you didn't follow any of this,
# go back, re-read, play with it, and ask for help. These are the foundations.


# ===   1.1.1  completing the database


# Next I'll add one more person, and create the other two tables:

x <- data.frame(id = autoincrement(philDB$person),
                name = "Kongzi",
                born = "551 BCE",
                died = "479 BCE",
                school = "Confucianism")
philDB$person <- rbind(philDB$person, x)

# a table of major works ...
philDB[["books"]] <- data.frame(id = 1:5,
                                title = c("Zhuangzi",
                                          "Analects",
                                          "Being and Time",
                                          "Daodejing",
                                          "On the Way to Language"),
                                published = c("300 BCE",
                                              "220 BCE",
                                              "1927",
                                              "530 BCE",
                                              "1959"))

# a "join" table that links works and their author ...
philDB[["works"]] <- data.frame(id = 1:5,
                                personID = c(3, 4, 2, 1, 2),
                                bookID = c(1, 2, 3, 4, 5))

str(philDB)


# ==   1.2  Querying the database  =============================================


# To retrieve data, we need to subset tables, possibly based on conditions we
# find in other tables. Sometimes we can simply get the information, e.g.
# all names ...
philDB$person$name
# ... or all book titles ...
philDB$books$title

# ... but sometimes we need to cross-reference information via join tables. Here
# is an example where we list authors and their works, sorted alphabetically by
# author:
(sel <- order(philDB$person$name))   # check out ?order and describe to
                                     # someone you know what it does, so that
                                     # you are sure you understand it. Its
                                     # indirection can be a bit tricky to
                                     # understand.
( pID <- philDB$person$id[sel] )
sel <- numeric()   # initialize the vector
for (ID in pID) {
  sel <- which(philDB$works$personID == ID)          # get all rows for which
                                                     # the condition is TRUE
  cat(sprintf("%s: ", philDB$person$name[ID]))       # output the person
  cat(sprintf("\"%s\"  ", philDB$books$title[sel]))  # output the book
  cat("\n")
}

# Examine the intermediate results and trace the logic until this is clear.


# ==   1.3  Task: Exercise  ====================================================

#    To make sure these principles are well understood, try the following
#    exercise:
#
#    Write code that adds another philosopher to the datamodel:
#       Immanuel Kant, (1724 - 1804), Enlightenment Philosophy.
#       Works: Critique of Pure Reason (1781), Critique of Judgement (1790)
#
#    Write code that lists the philosophical schools in
#    alphabetical order, and the books associated with them, also
#    alphabetically. Format your output like:
#    Confucianism
#       Analects - (220 BCE)
#    Daoism
#       Daodejing - (530 BCE)
#       ... etc.
#

# =    2  Implementing the protein datamodel  ==================================


# Working with the code above has probably illustrated the key concepts about
# curating data and storing it for analysis. In particular the join tables
# seem problematic - figuring out the correct IDs, it's easy to make
# mistakes.
#  - Data needs to be captured in a human-readable form so it can be verified
#      and validated;
#  - Some aspects of the database should _never_ be done by hand because
#      errors are easy to make and hard to see. That essentially includes
#      every operation that has to do with abstract, primary keys;
#  - Elementary operations we need to support are: adding data, selecting
#      data, modifying data and deleting data.

# We will therefore construct a protein database that we will use throughout
# the course as follows:
#  - For each table, we will keep the primary information in JSON files. JSON
#      files are easy to read, to edit if needed, and to modify.
#  - We will use simple scripts to read the JSON data and assemble it in
#      our database for further analysis.
#  - I have constructed initial files for yeast Mbp1 and nine other reference
#      species.
#  - I have written a small number of utility functions to read those files
#      and assemble them into a database. These are in .utilities.R which is
#      loaded on strtup of when this project.


# ==   2.1  JSON formatted source data  ========================================

# Have a look at the structure of the yeast Mbp1 protein data:
file.show("./data/MBP1_SACCE.json")

# - The whole thing is an array: [ ... ]. This is not necessary for a single
#     object, but we will have more objects in other files. And it's perfectly
#     legal to have an array with a single element.
# - The data is formatted as "key": "value" pairs inside an object { ... }.
#     This keeps the association between data items and their semantics
#     explicit.
# - All keys are strings and they are unique in the object.
# - Values are mostly single strings and integers ...
# - ... except for "sequence". That one is an array of strings. Why? This is to
#     make it easier to format and maintain the data. JSON does not allow line
#     breaks within strings, but the strings we copy/paste from Genbank or other
#     sources might have line breaks, sequence numbers etc. So we need to
#     sanitize the sequence at some point. But since we need to do that
#     anyway, it is easier to see the whole sequence if we store it in chunks.

# The .utilities.R script that get's loaded whenever you open this project
# has already made sure the "jsonlite" package exists on your computer. This
# package supports our work with .json formatted data.

if (! requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
# Package information:
#  library(help = jsonlite)       # basic information
#  browseVignettes("jsonlite")    # available vignettes
#  data(package = "jsonlite")     # available datasets

# try this out:
x <- jsonlite::fromJSON("./data/MBP1_SACCE.json")
str(x)

x$name
unlist(x$sequence)
sum(nchar(unlist(x$sequence))) # 833 amino acids


# ==   2.2  "Sanitizing" sequence data  ========================================


# Examine the dbSanitizeSequence() function:

dbSanitizeSequence

# Try:

dbSanitizeSequence(c("GAA", "ttc"))
dbSanitizeSequence("MsnQ00%0 I@#>YSary    S
                     G1 V2DV3Y>")
x <- "
        1 msnqiysary sgvdvyefih stgsimkrkk ddwvnathil kaanfakakr trilekevlk
       61 ethekvqggf gkyqgtwvpl niakqlaekf svydqlkplf dftqtdgsas pppapkhhha
      121 skvdrkkair sastsaimet krnnkkaeen qfqsskilgn ptaaprkrgr pvgstrgsrr
      ...
     " # this is copy/paste from Genbank
       # https://www.ncbi.nlm.nih.gov/protein/NP_010227.1/

x
dbSanitizeSequence(x)


# ==   2.3  Create a protein table for our data model  =========================

# ===   2.3.1  Initialize the database


# The function dbInit contains all the code to return a list of empty
# data frames for our data model.

dbInit

myDB <- dbInit()
str(myDB)


# ===   2.3.2  Add data


# fromJSON() returns a dataframe that we can readily process to add data
# to our table. Have a look at the function to add protein entries:

dbAddProtein

myDB <- dbAddProtein(myDB, jsonlite::fromJSON("./data/MBP1_SACCE.json"))
str(myDB)

# Lets check that the 833 amino acids of the yeast MBP1 sequence have
# safely arrived. Note the genral idiom we use here to retrieve the data:
# we define a boolean vector that satisfies a condition, then we subset
# a column with that vector.

sel <- myDB$protein$name == "MBP1_SACCE"
nchar(myDB$protein$sequence[sel])


# ==   2.4  Complete the database  =============================================


# Completing the database with Mbp1 data and data for 9 other "reference"
# species is more of the same. I have assembled the code in a script
# "./scripts/ABC-createRefDB.R" - open it, check it out, and then source it.
# It's really very simple, just reading some prepared files of data I have
# formatted with JSON, and assembling the data in our data model.
#
# The code is also very simple and in particular there is no checking for errors
# or inconsistencies. Have a look:

# Totally straightforward ...
dbAddTaxonomy
dbAddFeature

# Just slightly more complex, since we need to match the protein or feature
# name in the JSON file with its internal ID, and, when doing that confirm
# that it CAN be matched and that the match is UNIQUE
dbAddAnnotation

# Now: create the database
source("./scripts/ABC-createRefDB.R")

str(myDB)


# ===   2.4.1  Examples of navigating the database


# You can look at the contents of the tables in the usual way we access
# elements from lists and dataframes. Here are some examples:
myDB$protein
myDB$protein$RefSeqID
myDB$protein[,"name"]
myDB$taxonomy
myDB$taxonomy$species
biCode(myDB$taxonomy$species)

# Comparing two tables:
# Are all of the taxonomyIDs in the protein table present in the
# taxonomy table? We ought to check, because when we imported the
# data from JSON objects, we could have omitted or forgotten some. But we can
# check this with one simple expression. Unravel it and study its components.

all(myDB$protein$taxonomyID %in% myDB$taxonomy$ID)

# If this is not TRUE, you MUST fix the problem before continuing.

# Cross-referencing information:
# What is the species name of the protein whose name is "MBP1_COPCI"?

sel <- myDB$protein$name == "MBP1_COPCI"
x <- myDB$protein$taxonomyID[sel]

sel <- myDB$taxonomy$ID == x
myDB$taxonomy$species[sel]


# ==   2.5  Updating the database  =============================================


# Basic tasks for databases include retrieving data, selecting data, updating
# and deleting data. Here we will take a simple, pedestrian approach:
#
#     In case we need to modify any of the data, we modify it in the JSON file
#     save that, and recreate the database. The myDB database will only be
#     used for analysis.
#


# =    3  Add your own data  ===================================================


# You have been assigned a genome sequence fungus as "MYSPE", and your final task
# will be to find the protein in MYSPE that is most similar to yeast Mbp1, and
# to enter its information into the database.


# ==   3.1  Find a protein  ====================================================


# The BLAST algorithm will be properly introduced in a later learning unit -
# for now just use it in the following way:
#
# - Navigate to https://blast.ncbi.nlm.nih.gov/Blast.cgi and click on
#   Protein BLAST.
# - Enter NP_010227 into the "Query Sequence" field.
# - Choose "Reference proteins (refseq_protein)" as the "Database" in the
#   "Choose Search Set" section.
# - Paste your assigned MYSPE species name into the "Organism" field.
#
# - Click the "BLAST" button.

# You will probably get more than one result. If you get dozens of results or
# more, or if you get no results, something went wrong. Reconsider whether the
# problem was with your input, try something different, or ask for help.

# Otherwise, look for the top-hit in the "Descriptions" tab In some cases
# there will be more than one hit with nearly similar E-values. If this is the
# case for MYSPE, choose the one with the higher degree of similarity (more
# identities) with the N-terminus of the query - i.e. the Query sequence of
# the first ~ 100 amino acids.

# -  This is a crucial bit of coursework. Make sure it is well documented in
#    your journal.
# -  Follow the link to the protein data page, linked from "Accession".
# -  From there, in a separate tab, open the link to the taxonomy database page
#      for MYSPE which is linked from the "ORGANISM" record.


# ==   3.2  Put the information into JSON files  ===============================


# - Next make a copy of the file "./data/MBP1_SACCE.json" in the "data"
#     directory and give it a new name that corresponds to MYSPE - e.g. if
#     MYSPE is called "Crptycoccus neoformans", your file should be called
#     "MBP1_CRYNE.json"; in that case "MBP1_CRYNE" would also be the
#     "name" of your protein. Open the file in the RStudio editor and replace
#     all of the MBP1_SACCE data with the corresponding data of your protein.
#
#     Note: The UniProt ID may not be listed on the NCBI page. To retrieve
#     it, navigate to http://www.uniprot.org/id-mapping/ , paste your RefSeq ID
#     into the query field, make sure "RefSeqProtein" is selected for "From"
#     and "UniProtKB" is selected for "To", and click "Go". In case this does
#     not retrieve a single UniProt ID, contact me.
#
#     Save your .json file into your myScripts directory.
#
#     Confirm this step:
if (file.exists(sprintf("./myScripts/MBP1_%s.json", biCode(MYSPE)))) {
  cat("Excellent - all good to continue.\n")
} else {
  stop(sprintf(" The file \"./myScripts/MBP1_%s.json\" does not exist",
       biCode(MYSPE)))
}
#
#
# - Do a similar thing for the MYSPE taxonomy entry. Copy
#     "./data/refTaxonomy.json" and make a new file named "MYSPEtaxonomy.json".
#     Create a valid JSON file with only one single entry - that of MYSPE.
#
#     Confirm this step:
if (file.exists(sprintf("./myScripts/%staxonomy.json", biCode(MYSPE)))) {
  cat("Excellent - all good to continue.\n")
} else {
  stop(sprintf(" The file \"./myScripts/%staxonomy.json\" does not exist",
               biCode(MYSPE)))
}

# - Validate your two files online at https://jsonlint.com/


# ==   3.3  Create an R script to create your own database  ====================


# Next: to create your own database.
# - Make a new R script, call it "makeProteinDB.R"
# - enter the following expression as the first command:
#     source("./scripts/ABC-createRefDB.R")
# - than add the two commands that add your own protein and taxonomy data,
#     they should look like:
#
# myDB <- dbAddProtein(myDB,
#                      jsonlite::fromJSON("./myScripts/MBP1_<MYSPE>.json"))
# myDB <- dbAddTaxonomy(myDB,
#                       jsonlite::fromJSON("./myScripts/<MYSPE>taxonomy.json"))
#
#
# - save the .json file in the ./myScripts/ folder and source() it:
#
#     source("./myScripts/makeProteinDB.R")  # <<<- This command ...
#
# ... needs to be executed whenever you recreate the database. In particular,
# whenever you have added or modified data in any of the JSON files. Later you
# will add more information.

# Remember this principle. Don't rely on objects in memory - you might
# "break" them with a code experiment. But always have a script with
# which you can create what you need.


# ===   3.3.1  Check and validate

# Is your protein named according to the pattern "MBP1_MYSPE"? It should be.
# And does the taxonomy table contain the systematic name? It should be the same
# that you get when you type MYSPE into the console.

# Let's compute sequence lengths on the fly (with the function nchar() ), and
# open this with the table viewer function View()

View(cbind(myDB$protein[ , c("ID", "name", "RefSeqID")],
           length = nchar(myDB$protein$sequence)))

# Does your protein appear in the last row of this table? Where does your
# protein's length fall relative to the reference proteins? About the same? Much
# shorter? Much longer? If it is less then 500 amino acids long, I would suspect
# an error. Contact me for advice.

# Is that the right sequence? Is it the same as the one on the NCBI protein
# database page?
myDB$protein$sequence[nrow(myDB$protein)]

# If not, don't continue! Fix the problem first.
# Let me repeat: If this does not give you the right sequence of the MYSPE
#                Mbp1 homologue, DO NOT CONTINUE. Fix the problem.

# Is that the right taxonomy ID and binomial name for MYSPE?
# This question may be a bit non-trivial ... MYSPE is a species, but the
# recorded taxonomy ID may be a strain. We have a utility function,
# sameSpecies()  that normalizes organism name to the binomial species.
#
sel <- sameSpecies(myDB$taxonomy$species, MYSPE)
myDB$taxonomy[sel, ]

# If not, or if the result was "<0 rows> ... " then DO NOT CONTINUE.
# Fix the problem first.

# Does this give you the right refseq ID for MBP1_MYSPE?
sel <- myDB$protein$name == paste0("MBP1_", biCode(MYSPE))
myDB$protein$RefSeqID[sel]

# If not, or if the result was "<0 rows> ... " then DO NOT CONTINUE.
# Fix the problem first.


# ==   3.4  Record the results in your journal:  ===============================


# - In your journal, copy/paste the BLAST result headers from the
#     "Alignments" tab, to demonstrate that the data justifies your choice of
#     protein; you don't need to paste the whole alignment, just the header(s).
#     Note the relevant values separately: eValue, coverage, %ID etc. and link
#     to your protein's NCBI protein database page. (Note: in case there are
#     more than one high-scoring segments included for the SAME protein, you
#     need to show the results for all of its high-scoring segments.)
# - Copy and paste the contents of your two JSON files in your journal.
#     Make sure to format them in a fixed-width font, like Courier or Monaco
#     to keep them readable.
# - Here are two additional sanity checks.

# Has the MYSPE bicode been entered into the database?
biCode(myDB$taxonomy$species) %in% biCode(MYSPE)    # must be TRUE for the last
                                                    # element

sel <- sameSpecies(myDB$taxonomy$species, MYSPE)    # select where MYSPE species
                                                    # appears in the taxonomy
                                                    # table, and ...
myDB$protein$taxonomyID %in% myDB$taxonomy$ID[sel]  # confirm that a protein
                                                    # with the correct
                                                    # taxonomy ID exists in the
                                                    # table. Again, this must
                                                    # be TRUE for the last
                                                    # element.

# If these tests check out, your database is in a sane state and your creation
# script works as required.


# [END]
