# demoGsheet2R.R
#
# Demo to load data from a Google Sheet into R
#
# Demodata is here:
#
# 2022-09
# ==============================================================================
#
#

# A function to read Google sheet data via the REAST API:
readGsheet <- function (URL,
                        sheet = "") {

  # (1) Parse out the ID
  ID <- regmatches(URL, regexec("/d/([^/]+)/", URL))[[1]][2]

  # (2) make a retrieval URL for csv formatted data
  URL <- sprintf("https://docs.google.com/spreadsheets/d/%s%s%s%s",
                 ID,
                 "/gviz/tq?tqx=out:csv",
                 "&sheet=",
                 sheet)

  # (3) use the httr::GET() function to fetch the data.
  response <- httr::GET(URL)
  if (! httr::status_code(response) == 200) {
    stop(sprintf("Server status code was \"%s\".",
                 as.character(httr::status_code(response))))
  }

  # (4) extract the csv data from the response object
  x <-as.character(response)
  x <- strsplit(x, "\n")[[1]]
  x[1] <-      gsub("\\\"", "", x[1])

  # (5) read the csv text into an R data frame
  myDF <- read.csv(text = x)

  return(myDF)
}


# === DEMO   ===================================================================
if (FALSE) {

  myURL <- "https://docs.google.com/spreadsheets/d/1LGW-iV9egJlf4rEDBwfRumgmnmPr9PrJnkJn67mrS6Q/edit?usp=sharing"
  mySheet <- "Expenses"

  myExpenses <- readGsheet(myURL, mySheet)

  myCategories <- unique(myExpenses$Category)

  mySums <- numeric(length(myCategories))

  for (category in myCategories) {
    mySums[category] <- - sum(myExpenses$Amount[myExpenses$Category==category])
  }

  pie(mySums)

}


# [END]
