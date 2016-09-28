#' Test significant SNPs for bmass unit tests
#'
#' Test data created for checking bmass functionality
#' via unit tests, specifically a list of significant
#' SNPs
#'
#' @docType data
#'
#' @usage data(bmass_TestSigSNPs)
#'
#' @format A dataframe with 2 columns &  3 rows 
#' with the expected column headers `Chr` and `BP` 
#'
#' @keywords datasets
#'
#' @references Turchin MC and Stephens M. 2017. In preparation.
#'
#' @source Authors
#'
#' @examples
#' data(grav)
#' times <- attr(grav, "time")
#' phe <- grav$pheno
#' \donttest{iplotCurves(phe, times)}
"bmass_TestSigSNPs"
