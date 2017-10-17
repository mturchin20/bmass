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

#20170220 NOTE -- probably more official/better way to do below...check it out
library(ggplot2)
library(reshape2)

##Data1 <- read.table("../data/TestData1.txt", header=T)
##Data2 <- read.table("../data/TestData2.txt", header=T)
##SigSNPs <- read.table("../data/TestData1.GWASsnps.txt", header=T)
##source("PrepareData.R")
#library("devtools")
#devtools::load_all()
#devtools::test()
#data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
##bmass(c("Data1", "Data2"), GWASsnps=SigSNPs, NminThreshold = 2000, bmassSeedValue=NULL)
##bmassOutput1 <- bmass(c("Data1", "Data2"), GWASsnps=SigSNPs, NminThreshold = 2000, bmassSeedValue=NULL)
#bmassOutput1 <- bmass(c("bmass_TestData1", "bmass_TestData2"), GWASsnps=bmass_TestSigSNPs, NminThreshold = 2000, bmassSeedValue=NULL)
##ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N")
##DataList <- c("Data1", "Data2")

bmass <- function (DataSources, GWASsnps=NULL, MergedDataSources=NULL, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), GWASsnps_AnnotateWindow = 5e5, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), ProvidedPriors=NULL, UseFlatPriors=FALSE, GWASThreshFlag = 0, GWASThreshValue = 5e-8, NminThreshold = 0, PruneMarginalSNPs=TRUE, PruneMarginalSNPs_bpWindow=5e5, PrintMergedData=NULL, bmassSeedValue=NULL) {

#       print(DataSources)

        bmassOutput <- list()
	bmassOutput$MergedDataSources <- NULL
        if (!is.null(MergedDataSources)) {
		bmassOutput$MergedDataSources <- MergedDataSources
	}
	bmassOutput$ZScoresCorMatrix <- NULL
        bmassOutput$MarginalSNPs <- list()
	bmassOutput$Models <- NULL
	bmassOutput$ModelPriors <- NULL
        bmassOutput$PreviousSNPs <- list()
        bmassOutput$NewSNPs <- list()
        bmassOutput$GWASlogBFMinThreshold <- NULL
        #LogFile <- c()
	bmassOutput$LogFile <- c()

        bmassOutput$LogFile <- rbind(bmassOutput$LogFile, paste(format(Sys.time()), " -- beginning bmass.", sep=""))

        #20160823 CHECK_0: Prob -- list of Matthew functions specifically to double-check, go through, go over
        #       collapse
        #       em.priorprobs

        #Loading and checking data
        #~~~~~~

	if (!is.null(MergedDataSources)) {
		bmassOutput$LogFile <- rbind(bmassOutput$LogFile, paste(format(Sys.time()), " -- MergedDataSources was provided, skipping merging data step.", sep=""))
		
		#bmassOutput$LogFile <- CheckMergedDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, bmassOutput$MergedDataSources, ProvidedPriors, UseFlatPriors, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, bmassOutput$LogFile)
	}
	else {
		bmassOutput$LogFile <- CheckIndividualDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, bmassOutput$MergedDataSources, ProvidedPriors, UseFlatPriors, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, bmassOutput$LogFile)
		
		bmassOutput[c("MergedDataSources", "LogFile")] <- MergeDataSources(DataSources, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]
	}

	bmassOutput[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]

	bmassOutput[c("MergedDataSources", "MarginalSNPs", "ZScoresCorMatrix", "LogFile")] <- ProcessMergedAndAnnotatedDataSources(DataSources, bmassOutput$MergedDataSources, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, bmassOutput$LogFile)[c("MergedDataSources", "MarginalSNPs", "ZScoresCorMatrix", "LogFile")]

	if (is.null(PrintMergedData)) {
		bmassOutput$MergedDataSources <- NULL
	}

	bmassOutput[c("MarginalSNPs", "Models", "ModelPriors", "LogFile")] <- GetLogBFsFromData(DataSources, bmassOutput$MarginalSNPs, bmassOutput$ZScoresCorMatrix, SigmaAlphas, bmassOutput$LogFile)[c("MarginalSNPs", "Models", "ModelPriors", "LogFile")] 

	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "ModelPriors", "GWASlogBFMinThreshold", "LogFile")] <- DetermineAndApplyPriors(DataSources, bmassOutput$MarginalSNPs, GWASsnps, SigmaAlphas, bmassOutput$Models, bmassOutput$ModelPriors, ProvidedPriors, UseFlatPriors, GWASThreshFlag, GWASThreshValue, bmassSeedValue, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "ModelPriors", "GWASlogBFMinThreshold", "LogFile")]

	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")] <- FinalizeAndFormatResults(DataSources, bmassOutput$MarginalSNPs, bmassOutput$PreviousSNPs, GWASsnps, bmassOutput$GWASlogBFMinThreshold, SigmaAlphas, bmassOutput$Models, bmassOutput$ModelPriors, NminThreshold, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")]

#	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")] <- ExploreBestModelsUsingPosteriors(DataSources, bmassOutput$MarginalSNPs, bmassOutput$PreviousSNPs, bmassOutput$Models, bmassOutput$ModelPriors, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")]

##	bmassOutput$LogFile <- LogFile

	return(bmassOutput)

}


