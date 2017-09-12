# BIN-Storing_data.R
#
# Purpose: A Bioinformatics Course:
#              R code accompanying the BIN-Storing_data unit
#
# Version: 0.1
#
# Date:    2017  08  25
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 0.1    First code copied from BCH441_A03_makeYFOlist.R
#
# TODO:
#
#
# == HOW TO WORK WITH LEARNING UNIT FILES ======================================
#
# DO NOT SIMPLY  source()  THESE FILES!

# If there are portions you don't understand, use R's help system, Google for an
# answer, or ask your instructor. Don't continue if you don't understand what's
#  going on. That's not how it works ...
#
# ==============================================================================

# ==============================================================================
#        PART TWO: Implementing the (protein part of) the datamodel
# ==============================================================================

# Every entity in our datmodel can be implemented by a dataframe. To keep things
# organized, we create a list, that contains the enity tables. Here is a
# definition for a dataframe that holds the protein data for yeast Mbp1 as
# described in the data model schema:
#
# === Creating two tables ===
db <- list()      # we'll call our database "db" and start with an empty list

# Then add a dataframe with one row to the list
db$protein <- data.frame(ID = as.integer(1),
                         name = "Mbp1",
                         RefSeqID = "NP_010227",
                         UniProtID = "P39678",
                         taxonomy.ID = as.integer(4932),
                         sequence = "...",              # just a placeholder
                         stringsAsFactors = FALSE)

str(db)

# Next, we create the taxonomy table and add it into db
db$taxonomy <- data.frame(ID = as.integer(4932),
                          species = "Saccharomyces cerevisiae",
                          stringsAsFactors = FALSE)

# Let's add a second protein: We create a one-row dataframe with the data, then we rbind() it to the existing data frame:
myRow <- data.frame(ID = as.integer(2),
                    name = "Res2",
                    RefSeqID = "NP_593032",
                    UniProtID = "P41412",
                    taxonomy.ID = as.integer(4896),
                    sequence = "...",              # again, just a placeholder
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(4896),
                    species = "Schizosaccharomyces pombe",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

str(db)

# you can look at the contents of the tables in the usual way we would access
# elements from lists and dataframes. Here are some examples:
db$protein
db$protein$RefSeqID
db$protein[,"name"]
db$taxonomy
db$taxonomy$species
biCode(db$taxonomy$species)

# Here is an example to look up information in one table,
# based on a condition in another table:
# what is the species name for the protein
# whose name is "Mbp1"?

# First, get the taxonomy.ID for the Mbp1 protein. This is
# the key we need for the taxonomy table. We find it in a cell in the
# table: db$protein[<row>, <column>]
# <row> is that row for which the value in
# the "name" column is Mbp1:

db$protein$name == "Mbp1" # TRUE FALSE

# The <column> is called "taxonomy.ID". Simply
# insert these two expressions in the square
# brackets.

db$protein[db$protein$name == "Mbp1", "taxonomy.ID"]

# Assign the taxonomy.ID value ...
x <- db$protein[db$protein$name == "Mbp1", "taxonomy.ID"]

# ... and fetch the species_name value from db$taxonomy

db$taxonomy[db$taxonomy$ID == x, "species"]



# === Updating the database ====================================================

# Modifying a field

# Here is code that modifies the sequence field in the protein table of the
# database:
#
mySeq <- "
1 msnqiysary sgvdvyefih stgsimkrkk ddwvnathil kaanfakakr trilekevlk
61 ethekvqggf gkyqgtwvpl niakqlaekf svydqlkplf dftqtdgsas pppapkhhha
121 skvdrkkair sastsaimet krnnkkaeen qfqsskilgn ptaaprkrgr pvgstrgsrr
181 klgvnlqrsq sdmgfprpai pnssisttql psirstmgpq sptlgileee rhdsrqqqpq
241 qnnsaqfkei dledglssdv epsqqlqqvf nqntgfvpqq qssliqtqqt esmatsvsss
301 pslptspgdf adsnpfeerf pgggtspiis miprypvtsr pqtsdindkv nkylsklvdy
361 fisnemksnk slpqvllhpp phsapyidap idpelhtafh wacsmgnlpi aealyeagts
421 irstnsqgqt plmrsslfhn sytrrtfpri fqllhetvfd idsqsqtvih hivkrksttp
481 savyyldvvl skikdfspqy rielllntqd kngdtalhia skngdvvffn tlvkmgaltt
541 isnkegltan eimnqqyeqm miqngtnqhv nssntdlnih vntnnietkn dvnsmvimsp
601 vspsdyityp sqiatnisrn ipnvvnsmkq masiyndlhe qhdneikslq ktlksisktk
661 iqvslktlev lkesskdeng eaqtnddfei lsrlqeqntk klrkrliryk rlikqkleyr
721 qtvllnklie detqattnnt vekdnntler lelaqeltml qlqrknklss lvkkfednak
781 ihkyrriire gtemnieevd ssldvilqtl iannnknkga eqiitisnan sha
//
"


str(db$protein) # before
db$protein$sequence[db$protein$name == "Mbp1"] <- dbSanitizeSequence(mySeq)
str(db$protein) # after

# Analyze the expression ! Note how we specify an element of a vector (column)
# in a data frame in a list using a logical expression. And note how we assign
# the output (return value) of a function. As far as database coding goes this
# is pretty minimal - there is no error checking done at all. In particular: can
# we really guarantee that the name "Mbp1" is unique in the protein table? No!
# We never required it to be unique. This is a check we need to perform so
# frequently that we will encapsulate it in a function:

dbConfirmUnique

# try this:
dbConfirmUnique(c("TRUE", "FALSE"))
dbConfirmUnique(c(TRUE, FALSE))
dbConfirmUnique(c(TRUE, TRUE))
dbConfirmUnique(c(FALSE, FALSE))



#
#
# Here is the update to the sequence field of Res2 but using our
# confirmUnique() function


mySeq <- "
1 maprssavhv avysgvevye cfikgvsvmr rrrdswlnat qilkvadfdk pqrtrvlerq
61 vqigahekvq ggygkyqgtw vpfqrgvdla tkykvdgims pilsldideg kaiapkkkqt
121 kqkkpsvrgr rgrkpsslss stlhsvnekq pnssisptie ssmnkvnlpg aeeqvsatpl
181 paspnallsp ndntikpvee lgmleapldk yeeslldffl hpeegripsf lyspppdfqv
241 nsvidddght slhwacsmgh iemiklllra nadigvcnrl sqtplmrsvi ftnnydcqtf
301 gqvlellqst iyavdtngqs ifhhivqsts tpskvaaaky yldcilekli siqpfenvvr
361 lvnlqdsngd tslliaarng amdcvnslls ynanpsipnr qrrtaseyll eadkkphsll
421 qsnsnashsa fsfsgispai ispscsshaf vkaipsissk fsqlaeeyes qlrekeedli
481 ranrlkqdtl neisrtyqel tflqknnpty sqsmenlire aqetyqqlsk rlliwlearq
541 ifdlerslkp htslsisfps dflkkedgls lnndfkkpac nnvtnsdeye qlinkltslq
601 asrkkdtlyi rklyeelgid dtvnsyrrli amscginped lsleildave ealtrek
"

db$protein$sequence[dbConfirmUnique(db$protein$name == "Res2")] <- dbSanitizeSequence(mySeq)
str(db$protein)

# These expressions can get rather long, it's easier to read if we write:

select <- dbConfirmUnique(db$protein$name == "Res2")
value <- dbSanitizeSequence(mySeq)
db$protein$sequence[select] <- value

# ... and that's the code paradigm that we will adopt to update
# database fields (for now).

# Adding a record

# Adding a record is easy in principle, simply defining the values we get
# from the NCBI or EBI databases ... except for the ID field. That is a field
# we need to define internally, and once again, we'll use a small function
# to get this right.

dbAutoincrement

# After reading the code, predict the results
dbAutoincrement(1)      # Hm.
dbAutoincrement(1L)
dbAutoincrement(1:4)
dbAutoincrement(c(1, 2, 3, 4))
dbAutoincrement(c(1L, 2L, 3L, 4L))
dbAutoincrement(c(1, "4", 9))
dbAutoincrement(TRUE)

# Therefore, here is sample code to add one entry to the protein table.
#
mySeq <- "
1 msgdktifka tysgvpvyec iinnvavmrr rsddwlnatq ilkvvgldkp qrtrvlerei
61 qkgihekvqg gygkyqgtwi pldvaielae ryniqgllqp itsyvpsaad spppapkhti
121 stsnrskkii padpgalgrs rratsietes evigaapnnv segsmspsps dissssrtps
181 plpadrahpl hanhalagyn grdannhary adiildyfvt enttvpslli npppdfnpdm
241 sidddehtal hwacamgrir vvklllsaga difrvnsnqq talmratmfs nnydlrkfpe
301 lfellhrsil nidrndrtvf hhvvdlalsr gkphaaryym etminrlady gdqladilnf
361 qddegetplt maararskrl vrlllehgad pkirnkegkn aedyiieder frsspsrtgp
421 agielgadgl pvlptsslht seagqrtagr avtlmsnllh sladsydsei ntaekkltqa
481 hgllkqiqte iedsakvaea lhheaqgvde erkrvdslql alkhainkra rddlerrwse
541 gkqaikrarl qaglepgals tsnatnapat gdqkskddak sliealpagt nvktaiaelr
601 kqlsqvqank telvdkfvar areqgtgrtm aayrrliaag cggiapdevd avvgvlcell
661 qeshtgarag aggerddrar dvammlkgag aaalaanaga p
"

myRow <- data.frame(ID = dbAutoincrement(db$protein$ID),
                    name = "UMAG_1122",
                    RefSeqID = "XP_011392621",
                    UniProtID = "A0A0D1DP35",
                    taxonomy.ID = as.integer(5270),
                    sequence = dbSanitizeSequence(mySeq),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(5270),
                    species = "Ustilago maydis",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

# Here is a bit of code template for the same ... it won't execute of course but
# you can copy it into your own code script file to modify when you need to add
# your own protein.

# ===== begin code template =====================================
mySeq <- "
RAW SEQUENCE FROM NCBI RECORD
"

myRow <- data.frame(ID = dbAutoincrement(db$protein$ID),
                    name = "NAME",
                    RefSeqID = "REFSEQ",
                    UniProtID = "UNIPROT",
                    taxonomy.ID = as.integer(TAX_ID),
                    sequence = dbSanitizeSequence(mySeq),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(TAX_ID),
                    species = "BINOMIAL NAME",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

# ===== end code template ===========================================




# deleting a record

# This is simple code without error checking and of course we can make mistakes.
# Often we can just overwrite the offending field with correct data. But
# sometimes it will be easier (and more robust) to delete the erroneous entry
# and add the correct one. For example, if your code is in a script, and you
# realize the entry had an error, I would not "patch" the error in the script
# but delete the entry on the console command line, correct the script and
# execute the correct block. That way we fix the mistake at its source.
#
# Removing a row from a datframe is trivial: just overwrite the dataframe with
# a selection statement in which the unique selection of the offending row is
# inverted:

# create an erroneous entry
myRow <- data.frame(ID = dbAutoincrement(db$protein$ID),
                    name = "nonesuch",
                    RefSeqID = "NP_000000",
                    UniProtID = "A123456",
                    taxonomy.ID = as.integer(999999),
                    sequence = dbSanitizeSequence("utter nonsense"),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

# Make a logical vector that identifies it
select <- dbConfirmUnique(db$protein$name == "nonesuch")
select     # the selection
!select    # its logical inversion

str(db)    # before
db$protein <- db$protein[ ! select, ]  # overwrite the table with a copy
# without the selected record
str(db)    # after

# Note: if you delete records "by hand" you need to be careful that you do
# not remove keys that are used as foreign keys in another table - if there
# are such dependecies, you need to update the other table(s) too. "Real"
# database systems include such dependencies in the creation instructions
# of the table schema: "on delete cascade ..."

#
# ==============================================================================



# ==============================================================================
#        PART ONE: Database maintenance
# ==============================================================================

# === Merging databases ========================================================

# From time to time, we will want to update the systems database - either
# because the schema has changed, e.g. we may include additional tables, or
# because I have added reference data: sequences, annotations and the like. The
# goal here is to make this as simple as possible, keeping all information you
# may have entered intact, and avoiding to duplicate entries. This is not quite
# as easy as one would hope. First of all, since we need to merge records from
# different databases - your own, and the reference - we can't assume that the
# primary keys are unique. You may have introduced the same proteinID that I
# used later for the reference database. Thus we need to update keys. But if we
# update primary keys, we need to make sure that every record that references a
# key as a foreign key is updated as well. And ensuring database integrity for
# "cascading updates" takes a _lot_ of care to guarantee correctness. You don't
# want to read through that code. I don't want to write it. And we should not
# end up pretending to build a real and proper database engine here!
#
# We'll do something simple instead. The sanest way is to create
# primary keys with "namespaces", ie. the keys will reflect the source of the
# information. This way, even though the databases are merged, the information
# sources can be identified. Still, it can be accessed through a common
# mechanism. A database merge with separate namespaces can simply proceed as
# follows:

#
#    for each table in the database
#        merge the table from both source databases
#        throw out any row with a duplicate ID
#
# In this way, if we merge your working database with the reference database
# all reference database records will be updated, and all your working records
# will remain.

# And while we are adding semnatics to keys, we should also identify in which
# table they are being used. You may remember that I told you: keys should not
# try to express semantics? Here we have a slightly different situation: since
# we are not actually constructing an underlying database engine, we will often
# encounter keys that we use "manually". And in that case it is very helpful to
# make the key's domain explicit.

# The code I wrote to implement this was loaded when you typed init().

# Let's have a brief look at the new dbAutoincrement() function. You'll see the
# code as usual when you execute the function name without its parentheses:

dbAutoincrement

# The default namespace is "my" - that's the namespace you use. I will use
# "ref" for all entries of the reference database. Thus you normally won't
# need to enter the namespace identifier - only I do that.

refDB$protein$ID
dbAutoincrement(refDB$protein$ID)
dbAutoincrement(refDB$protein$ID, ns = "ref")

# Once the namespaces are thus kept well apart, we won't overwrite each other's
# work. Merging two table becomes a breeze:

tableUnion

# Here is an example
this <- data.frame(ID = "ref_1", type = "one crow", stringsAsFactors = FALSE)
this <- rbind(this, data.frame(ID = dbAutoincrement(this$ID, ns="ref"),
                               type = "for",
                               stringsAsFactors = FALSE))
this <- rbind(this, data.frame(ID = dbAutoincrement(this$ID, ns="ref"),
                               type = "sorrow",
                               stringsAsFactors = FALSE))
that <- data.frame(ID = "my_1", type = "two crows", stringsAsFactors = FALSE)
that <- rbind(that, this[2:3, ])  # repeat some rows
that <- rbind(that, data.frame(ID = dbAutoincrement(that$ID),
                               type = "for ",
                               stringsAsFactors = FALSE))
that <- rbind(that, data.frame(ID = dbAutoincrement(that$ID),
                               type = "joy ...",
                               stringsAsFactors = FALSE))

that
this

rhyme <- tableUnion(this, that)
rhyme$type

# Finally, we need to do this for all tables in the datamodel. Well, almost all:
# we introduce a field $version, through which we can ensure that the datamodels
# have the same schema. If they would have different schemas the merging would
# break. But we don't merge this single field. (Incidentally: a schema version
# is not "data", rather we call this "metadata": information _about_ the data.)
# So we exclude the "version" from the names of elements to merge.
#
# But how to iterate over all tables in a list by name? You might be tempted to
# try something like
#
# n <- names(db)
# myDB$n[1], myDB$n[2] etc.
# ... but NO! That's not how the $ operator works.
#
# The R community keeps particularly poignant comments from the R-help mailing
# list in a package called "fortunes", and fortune(312) reads:
#
# " The problem here is that the $ notation is a magical shortcut and like any
#   other magic if used incorrectly is likely to do the programmatic equivalent
#   of turning yourself into a toad.
#
# -- Greg Snow (in response to a user that wanted
#    to access a column whose name is stored in y
#    via x$y rather than x[[y]])
#    R-help (February 2012)

# So in order to avoid turning into toads, we create a vector "tables", iterate
# over "table" elements from "tables" and use them as ref[[table]] ... etc.

# === Putting data into the new schema =========================================

# Entering the YFO data into the new schema takes almost exactly
# the same syntax as the code you wrote for the last assignment.
#

# === TASK: ===
# Initialize myDB with the YFO homologue of yeast Mbp1

# First: initialize "myDB", this is the name of the R object that
#        will contain data you collect.

myDB <- dbInit()   # This creates an empty database with the latest schema

# Next: Add your protein and the YFO to the database. Copy the code-template
#       below to your myCode.R file, edit it to replace the placeholder
#       items with your data, and execute the code. Then continue below
#       the code template.


# ===== begin code template: add a protein and an organism to the database =====

# == edit placeholder items!
myBinomial <- "<BINOMIAL NAME>"    # Name, NOT biCode()
myTaxonomyId <- as.integer(<TAX_ID>)
myProteinName <- "<PROTEIN NAME>"  # Name your protein "MBP1_<YFO>" (where <YFO>
# is a placeholder for the biCode() of YFO)
myProteinRefSeqID <- "<REFSEQID>"
myProteinUniProtID <- "<UNIPROTID>"
mySeq <- "
<SEQUENCE>
"

# == create the protein entry
proteinRow <- data.frame(ID = dbAutoincrement(myDB$protein$ID),
                         name = myProteinName,
                         RefSeqID = myProteinRefSeqID,
                         UniProtID = myProteinUniProtID,
                         taxonomy.ID = myTaxonomyId,
                         sequence = dbSanitizeSequence(mySeq),
                         stringsAsFactors = FALSE)
myDB$protein <- rbind(myDB$protein, proteinRow)

# == create the taxonomy entry
taxonomyRow <- data.frame(ID = myTaxonomyId,
                          species = myBinomial,
                          stringsAsFactors = FALSE)
myDB$taxonomy <- rbind(myDB$taxonomy, taxonomyRow)
# ===== end code template ===========================================

# ... continue here.

# myDB now contains one record each in two tables. The remaining tables exist
# but they are empty. Check that all the values are correct: just execute
myDB

# Now let's merge myDB with the data from refDB. refDB should already have been
# loaded from .utilities.R ; you can also explore the original script with which
# refDB was created, for your reference: it is create_refDB.R  The database
# refDB is the default argument for dbMerge(), so you don't need to
# specify it. By assigning the result of dbMerge() back to myDB we overwrite the
# previous version.

myDB <- dbMerge(myDB)

str(myDB)

# check the protein table
View(myDB$protein[ , c("ID", "name", "RefSeqID")])

# Is your protein named according to the pattern "MBP1_<YFO>"? It should be.
# And does the taxonomy table contain the binomial name? It should be the same
# that you get when you type YFO into the console.

# Let's compute sequence lengths on the fly (with the function nchar() ) and
# add them to our view. Then we'll open this with the table viewer function
# View()

View(cbind(myDB$protein[ , c("ID", "name", "RefSeqID")],
           length = nchar(myDB$protein$sequence)))

# Where does your protein's length fall relative to the reference proteins?
# About the same? Much shorter? Much longer?

# Is that the right sequence?
sel <- myDB$protein$ID == "my_pro_1"
myDB$protein$sequence[sel]

# If not, don't continue! Fix the problem first.
# Let me repeat: If this does not give you the right sequence of the YFO
#                Mbp1 homologue, DO NOT CONTINUE. Fix the problem.

# Is that the right taxonomy ID and binomial name for YFO?
sel <- myDB$taxonomy$species == YFO
myDB$taxonomy[sel, ]

# If not, or if the result was "<0 rows> ... " then DO NOT CONTINUE.
# Fix the problem first.

# Does this give you the right protein ID for MBP1_<YFO>?
sel <- dbConfirmUnique(myDB$protein$name == paste("MBP1_", biCode(YFO), sep = ""))
myDB$protein$ID[sel]

# If not, or if the result was "<0 rows> ... " then DO NOT CONTINUE.
# Fix the problem first.

#
#
# === Saving and loading data ==================================================
#
# Once you have confirmed that myDB has been successfully created and updated
# and is in a good state, you should make a backup copy. There are many ways to
# save data to a file on disk and read it back in. One of the most convenient is
# the function pair save() and load().
#
# save() saves an R object to a file. Its signature is

#    save(<object-name/s>, file = <file-name>)

# The object, or objects that are saved like this are identically recreated when
# load() is executed. The object name is not changed - you don't assign the
# result of load(). You use load() for its "side-effect": re-creating the saved
# object, not for using the return-value. Therefore the signature for load() is
# simply

#     load(<file-name>)

# All R objects in <file-name> are created by load(), if they don't yet exist.
# If they already exist, the will be overwritten. The only effect you see is
# that the object appears in the Environment pane of RStudio (if it wasn't there
# already), or you may notice that its content has changed if it did exist.

# == TASK: ==
# Save myDB so you can recreate it easily when you next open RStudio.

save(myDB, file = "myDB.01.RData")  # Note that I give this a version number!

# Let's confirm:
rm(myDB)
nrow(myDB$protein)  # must give an error
load("myDB.01.RData")
nrow(myDB$protein)  # must be 11. If not, don't continue etc.

===New Database ===

  Here is some sample code to work with the new database, enter new protein data for YFO, save it and load it again when needed.


<source lang="R">
  # You don't need to load the reference database refDB. If
  # everything is set up correctly, it gets loaded at startup.
  # (Just so you know: you can turn off that behaviour if you
  # ever should want to...)


  # First you need to load the newest version of dbUtilities.R

  updateDButilities("7bb32ab3d0861ad81bdcb9294f0f6a737b820bf9")

# If you get an error:
#    Error: could not find function "updateDButilities"
# ... then it seems you didn't do the previous update.

# Try getting the update with the new key but the previous function:
# updateDbUtilities()
#
# If that function is not found either, confirm that your ~/.Rprofile
# actually loads dbUtilites.R from your project directory.

# As a desperate last resort, you could uncomment
# the following piece of code and run the update
# without verification...
#
# URL <- "http://steipe.biochemistry.utoronto.ca/abc/images/f/f9/DbUtilities.R"
# download.file(URL, paste(PROJECTDIR, "dbUtilities.R", sep="")), method="auto")
# source(paste(PROJECTDIR, "dbUtilities.R", sep=""))
#
# But be cautious: there is no verification. You yourself need
# to satisfy yourself that this "file from the internet" is what
# it should be, before source()'ing it...


# After the file has been source()'d,  refDB exists.
ls(refDB)


# check the contents of refDB:
refDB$protein$name
refDB$taxonomy


# list refSeqIDs for saccharomyces cerevisiae genes.
refDB$protein[refDB$protein$taxID == 559292, "refSeqID"]


# To add some genes from YFO, I proceed as follows.
# Obviously, you need to adapt this to your YFO
# and the sequences in YFO that you have found
# with your PSI-BLAST search.

# Let's assume my YFO is the fly agaric (amanita muscaria)
# and its APSES domain proteins have the following IDs
# (these are not refSeq btw. and thus unlikely
# to be found in UniProt) ...
# KIL68212
# KIL69256
# KIL65817
#


# First, I create a copy of the database with a name that
# I will recognize to be associated with my YFO.
amamuDB <- refDB


# Then I fetch my protein data ...
tmp1 <- fetchProteinData("KIL68212")
tmp2 <- fetchProteinData("KIL69256")
tmp3 <- fetchProteinData("KIL65817")


# ... and if I am satisfied that it contains what I
# want, I add it to the database.
amamuDB <- addToDB(amamuDB, tmp1)
amamuDB <- addToDB(amamuDB, tmp2)
amamuDB <- addToDB(amamuDB, tmp3)


# Then I make a local backup copy. Note the filename and
# version number  :-)
save(amamuDB, file="amamuDB.01.RData")


# Now I can explore my new database ...
amamuDB$protein[amamuDB$protein$taxID == 946122, "refSeqID"]


# ... but if anything goes wrong, for example
# if I make a mistake in checking which
# rows contain taxID 946122 ...
amamuDB$protein$taxID = 946122

# Ooops ... what did I just do wrong?
#       ... wnat happened instead?

amamuDB$protein$taxID


# ... I can simply recover from my backup copy:
load("amamuDB.01.RData")
amamuDB$protein$taxID






# [END]
