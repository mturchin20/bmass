#' Test data 1 for bmass unit tests
#'
#' Test data created for checking bmass functionality
#' via unit tests.
#'
#' @docType data
#'
#' @usage data(bmass_TestData1)
#'
#' @format A dataframe with 7 columns &  11 rows 
#' with the expected column headers for bmass
#'
#' @keywords datasets
#'
#' @references Turchin MC and Stephens M. 2017. In preparation.
#'
#' @source Authors
#'
#' @references Moore et al. (2013) Genetics 195:1077-1086
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/23979570}{PubMed})
#'
#' @source \href{http://qtlarchive.org/db/q?pg=projdetails&proj=moore_2013b}{QTL Archive}
#'
#' @examples
#' data(grav)
#' times <- attr(grav, "time")
#' phe <- grav$pheno
#' \donttest{iplotCurves(phe, times)}
"bmass_TestData1"
