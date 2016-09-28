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

##Data1 <- read.table("../data/TestData1.txt", header=T)
##Data2 <- read.table("../data/TestData2.txt", header=T)
##SigSNPs <- read.table("../data/TestData1.GWASsnps.txt", header=T)
##source("PrepareData.R")
#library("devtools")
#devtools::load_all()
#data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
##bmass(c("Data1", "Data2"), GWASsnps=SigSNPs, NminThreshold = 2000, bmassSeedValue=NULL)
#bmassOutput1 <- bmass(c("Data1", "Data2"), GWASsnps=SigSNPs, NminThreshold = 2000, bmassSeedValue=NULL)
##ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N")
##DataList <- c("Data1", "Data2")


#bmass <- function(ExpectedColumnNames, DataFileNames, DataFileLocations, OutputFileBase) {
#
#	PrepareData(ExpectedColumnNames, DataFileNames, DataFileLocations, OutputFileBase)
#
#}

bmass <- function (DataSources, GWASsnps=NULL, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), MergedDataSources=NULL, ProvidedPriors=NULL, UseFlatPriors=FALSE, PruneMarginalHits=TRUE, PruneMarginalHits_bpWindow=5e5, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, NminThreshold = 0, bmassSeedValue=NULL) {

        print(DataSources)

        LogFile <- c()
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

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- beginning bmass.", sep=""))

        #20160823 CHECK_0: Prob -- list of Matthew functions specifically to double-check, go through, go over
        #       collapse
        #       em.priorprobs

        #Loading and checking data
        #~~~~~~

	LogFile <- CheckIndividualDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile)

#	<- MergeDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile)

#	<- AnnotateMergedDataWithGWASSNPs(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile)

#	<- ProcessMergedAndAnnotatedDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile)

#	<- GetLogBFsFromData(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile)

#	<- FinalizeAndFormatResults(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile)

	bmassOutput$LogFile <- LogFile

	return(bmassOutput)

}




#if (FALSE) {
#~~~
#TimeTag 20160927 -- Starting to develop and connect pieces of new, subdivided .R code files


#}

