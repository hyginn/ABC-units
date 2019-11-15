# ABC-createRefDB.R
#
# Create a reference protein database for Mbp1-like proteins
#
# Boris Steipe for ABC learning units
#
# For the species, see:
# http://steipe.biochemistry.utoronto.ca/abc/index.php/Reference_species_for_fungi
#
# For the data model, see
# https://docs.google.com/drawings/d/1uupNvz18_FYFwyyVPebTM0CUxcJCPDQuxuIJGpjWQWg
# For the schema, see dbInit() in ./scripts/ABC-dbUtilities.R
#
# ==============================================================================


myDB <- dbInit()

myDB <- dbAddProtein(myDB, jsonlite::fromJSON("./data/MBP1_SACCE.json"))
myDB <- dbAddProtein(myDB, jsonlite::fromJSON("./data/refMBP1Proteins.json"))
myDB <- dbAddProtein(myDB, jsonlite::fromJSON("./data/refAPSES_PSI-BLAST.json"))

myDB <- dbAddTaxonomy(myDB, jsonlite::fromJSON("./data/refTaxonomy.json"))

myDB <- dbAddFeature(myDB, jsonlite::fromJSON("./data/refFeatures.json"))

myDB <- dbAddAnnotation( myDB, jsonlite::fromJSON("./data/refAnnotations.json"))


# [END]
