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

myDB <- dbAddProtein(    myDB, fromJSON("./data/MBP1_SACCE.json"))
myDB <- dbAddProtein(    myDB, fromJSON("./data/refProteins.json"))
myDB <- dbAddTaxonomy(   myDB, fromJSON("./data/refTaxonomy.json"))
myDB <- dbAddFeature(    myDB, fromJSON("./data/refFeatures.json"))
myDB <- dbAddAnnotation( myDB, fromJSON("./data/refAnnotations.json"))


# [END]
