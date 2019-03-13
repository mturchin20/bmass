#' @name bmass_TestData1
#'
#' @title bmass Example Dataset 1
#'
#' @docType data
#' 
#' @description A manually created sample dataset for use in vignettes
#' and unit tests.
#'
#' @format A data frame with 11 rows and 9 variables:
#' 
#' \describe{
#'   \item{Chr}{chromosome}
#'   \item{BP}{basepair position}
#'   \item{Marker}{rsID# or other identifier}
#'   \item{MAF}{Minor Allele Frequency}
#'   \item{A1}{reference allele}
#'   \item{A2}{alternative allele}
#'   \item{Direction}{direction of association
#'	effect size, + or -}
#'   \item{pValue}{p-Value of GWAS association}
#'   \item{N}{sample size}
#' }
#' 
#' @source Manually created
#'
#' @keywords data
#' 
NULL

#' @name bmass_TestData2
#'
#' @title bmass Example Dataset 2
#'
#' @docType data
#'
#' @description A manually created sample dataset for use in vignettes
#' and unit tests.
#'
#' @format A data frame with 11 rows and 9 variables:
#' 
#' \describe{
#'   \item{Chr}{chromosome}
#'   \item{BP}{basepair position}
#'   \item{Marker}{rsID# or other identifier}
#'   \item{MAF}{Minor Allele Frequency}
#'   \item{A1}{reference allele}
#'   \item{A2}{alternative allele}
#'   \item{Direction}{direction of association
#'	effect size, + or -}
#'   \item{pValue}{p-Value of GWAS association}
#'   \item{N}{sample size}
#' }
#'
#' @keywords data
#' 
#' @source Manually created
#'
NULL

#' @name bmass_TestSigSNPs
#' 
#' @title bmass Example GWAS SNPs
#'
#' @docType data
#' 
#' @description A manually created list of GWAS significant SNPs to be
#' used in conjunction with 'bmass_TestData1' and 'bmass_TestData2'.
#'
#' @format A data frame with 2 rows and 2 variables:
#' 
#' \describe{
#'   \item{Chr}{chromosome}
#'   \item{BP}{basepair position}
#' }
#'
#' @keywords data
#' 
#' @source Manually created
#' 
NULL

#' @name GlobalLipids2013.GWASsnps
#'
#' @title GlobalLipids2013 GWAS SNPs 
#'
#' @docType data
#' 
#' @description A list of the univariate GWAS significant SNPs from
#' the GlobalLipids2013 dataset to be used in the main bmass vignette.
#'
#' @format A data frame with 157 rows and 2 variables:
#' 
#' \describe{
#'   \item{Chr}{chromosome}
#'   \item{BP}{basepair position}
#' }
#' 
#' @source Supplementary Tables 2 and 3 from
#' \url{https://www.nature.com/articles/ng.2797}.
#' 
#' @keywords data
#'
NULL
