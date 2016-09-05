#' Checking and setting up logfile interface.
#'
#' Description
#'
#' @param x Something.
#' @param y Something2.
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)





#bmass <- function(ExpectedColumnNames, DataFileNames, DataFileLocations, OutputFileBase) {
#
#	PrepareData(ExpectedColumnNames, DataFileNames, DataFileLocations, OutputFileBase)
#
#}

bmass <- function (DataSources, GWASsnps=NULL, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), MergedDataSources=NULL, ProvidedPriors=NULL, UseFlatPriors=FALSE, PruneMarginalHits=TRUE, PruneMarginalHits_bpWindow=5e5, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, NminThreshold = 0, bmassSeedValue=NULL) {

        print(DataSources)

        LogFile1 <- c()
        bmassOutput <- list()
        bmassOutput$ModelPriors <- NULL
        bmassOutput$MarginalSNPs <- list()
        bmassOutput$MarginalSNPs$SNPs <- NULL
        bmassOutput$MarginalSNPs$logBFs <- NULL
        bmassOutput$MarginalSNPs$Posteriors <- NULL
        bmassOutput$NewSNPs <- list()
        bmassOutput$NewSNPs$SNPs <- NULL
        bmassOutput$NewSNPs$logBFs <- NULL
        bmassOutput$NewSNPs$Posteriors <- NULL
        bmassOutput$GWASlogBFMinThreshold <- NULL
        bmassOutput$PreviousSNPs <- list()
        bmassOutput$PreviousSNPs$SNPs <- NULL
        bmassOutput$PreviousSNPs$logBFs <- NULL
        bmassOutput$PreviousSNPs$Posteriors <- NULL

        LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning bmass.", sep=""))

        #20160823 CHECK_0: Prob -- list of Matthew functions specifically to double-check, go through, go over
        #       collapse
        #       em.priorprobs

        #Loading and checking data
        #~~~~~~

	LogFile1 <- DataChecks(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile1)

}










