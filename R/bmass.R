#' @title Bayesian multivariate analysis of summary statistics 
#' (\code{bmass})
#'
#' @description Run \code{bmass} on a set of phenotypes that each have
#' univariate GWAS statistics on the same set of SNPs
#' 
#' @param DataSources A string indicating the variable names of the
#' input datafiles and phenotypes. No default value.
#' 
#' @param GWASsnps A data.table containing rows of SNPs that were
#' univariate genome-wide significant in the phenotypes being used for
#' analysis; \code{GWASsnps} input file should have two columns, one 
#' for chromosome and another for basepair position (with column headers
#' of \code{Chr} and \code{BP}). No default value.
#' 
#' @param GWASThreshFlag A logical \code{TRUE}/\code{FALSE} flag that 
#' indicates whether to threshold input \code{GWASsnps} list by a 
#' univariate GWAS p-value or not (eg the input \code{GWASsnps} list 
#' contains variants that are significant from discovery + replication 
#' data, but the input summary statistics are just from the discovery 
#' cohort). Default is \code{TRUE}.
#' 
#' @param GWASThreshValue A numerical value indicating the univariate
#' p-value threshold to use in conjunction with the \code{GWASThreshFlag}.
#' Default is \code{5e-8}.
#' 
#' @param SNPMarginalUnivariateThreshold A numerical value indicating
#' the univariate p-value threshold to use when collecting marginally
#' significant SNPs for final \code{bmass} analysis. Default is 
#' \code{1e-6}.
#' 
#' @param SNPMarginalMultivariateThreshold A numerical value
#' indicating the basic multivariate p-value threshold to use when
#' collecting marginally significant SNPs for final \code{bmass} 
#' analysis. Default is \code{1e-6}.
#' 
#' @param NminThreshold A numerical value that indicates a sample size
#' threshold to use where SNPs below which are removed. Default is
#' \code{0}. 
#' 
#' @param PrintMergedData A logical \code{TRUE}/\code{FALSE} flag that 
#' indicates whether the intermediary 'merged datafile' should be included 
#' in the final \code{bmass} output; this file combines all the phenotypes 
#' for every SNP provided just prior to thresholding for marginally
#' significant SNPs. Default is \code{FALSE}.
#' 
#' @param PrintProgress A logical \code{TRUE}/\code{FALSE} flag that 
#' indicates whether progress statements should be printed to 
#' \code{stderr} during the course of running \code{bmass} or not.
#' Default is \code{FALSE}.
#' 
#' @param ... Additional optional arguments.
#' 
#' @return A list containing model, SNP, and posterior information for
#' both the previously significant univariate SNPs (\code{PreviousSNPs})
#' and the newly significant multivariate SNPs (\code{NewSNPs}). For a 
#' full breakdown of the \code{bmass} output list structure, please see 
#' the associated vignettes.
#'
#' @examples
#' Phenotypes <- c("bmass_SimulatedData1", "bmass_SimulatedData2")
#' bmassOutput <- bmass(Phenotypes, bmass_SimulatedSigSNPs)
#' summary(bmassOutput)
#' bmassOutput$NewSNPs$SNPs
#'
#' @section Other Examples:
#' \code{bmass(c("HDL","LDL","TG","TC"), GWASsnps, NminThreshold = 50000)}
#' \code{bmass(c("HDL","LDL","TG","TC"), GWASsnps, GWASThreshValue = 1e-8,
#'   NminThreshold = 50000, PrintProgress = TRUE)} 
#' \code{bmass(c("HDL", "LDL", "TG", "TC"), GWASsnps, GWASThreshFlag = FALSE,
#'   SNPMarginalUnivariateThreshold = 1e-4,
#'   SNPMarginalMultivariateThreshold = 1e-4,
#'   PrintMergedData = TRUE)} 
#' \code{bmassOutput <- bmass(c("HDL","LDL","TG","TC"),
#'   GWASsnps, NminThreshold = 50000)} 
#'
#' @export
#' 
bmass <- function (DataSources, GWASsnps = NULL, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, GWASThreshFlag = TRUE, GWASThreshValue = 5e-8, NminThreshold = 0, PrintMergedData = FALSE, PrintProgress = FALSE, ...) {
	return(bmassMain(DataSources, GWASsnps, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, GWASThreshFlag, GWASThreshValue, NminThreshold, PrintMergedData, PrintProgress, ...))
}

bmassMain <- function (DataSources, GWASsnps, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, GWASThreshFlag, GWASThreshValue, NminThreshold, PrintMergedData, PrintProgress, MergedDataSources = NULL, ZScoresCorMatrix = NULL, ExpectedColumnNames = c("Chr", "BP", "Marker", "MAF", "A1", "Direction", "pValue", "N"), GWASsnps_AnnotateWindow = 5e5, SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), ProvidedPriors = NULL, UseFlatPriors = NULL, PruneMarginalSNPs = TRUE, PruneMarginalSNPs_bpWindow = 5e5, bmassSeedValue = 1) {

        bmassOutput <- list()
	bmassOutput$MergedDataSources <- NULL
        if (!is.null(MergedDataSources)) {
		bmassOutput$MergedDataSources <- MergedDataSources
	}
	bmassOutput$ZScoresCorMatrix <- NULL
        if (!is.null(ZScoresCorMatrix)) {
		bmassOutput$ZScoresCorMatrix <- ZScoresCorMatrix
	}
	bmassOutput$MarginalSNPs <- list()
	bmassOutput$Models <- NULL
	bmassOutput$ModelPriors <- NULL
        bmassOutput$PreviousSNPs <- list()
        bmassOutput$NewSNPs <- list()
        bmassOutput$GWASlogBFMinThreshold <- NULL
	bmassOutput$LogFile <- c()

        bmassOutput$LogFile <- rbind(bmassOutput$LogFile, paste(format(Sys.time()), " -- beginning bmass.", sep=""))
	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- beginning bmass.", sep=""), stderr())
	}

        #Loading and checking data
        #~~~~~~

	if (!is.null(MergedDataSources)) {
		bmassOutput$LogFile <- rbind(bmassOutput$LogFile, paste(format(Sys.time()), " -- MergedDataSources was provided, skipping merging data step.", sep=""))
		if (PrintProgress == TRUE) {
        		write(paste(format(Sys.time()), " -- MergedDataSources was provided, skipping merging data step.", sep=""), stderr())
		}
	} else {
		if (PrintProgress == TRUE) {
        		write(paste(format(Sys.time()), " -- Checking individual datasource files and merging datasets.", sep=""), stderr())
		}
		bmassOutput$LogFile <- CheckIndividualDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, bmassOutput$MergedDataSources, ProvidedPriors, UseFlatPriors, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, bmassOutput$LogFile)
		
		bmassOutput[c("MergedDataSources", "LogFile")] <- MergeDataSources(DataSources, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]
	}

	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Annotating merged datasources with GWAS SNPs (if provided) then processing finalized merged dataset.", sep=""), stderr())
	}
	bmassOutput[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]

	bmassOutput[c("MergedDataSources", "MarginalSNPs", "ZScoresCorMatrix", "LogFile")] <- ProcessMergedAndAnnotatedDataSources(DataSources, bmassOutput$MergedDataSources, bmassOutput$ZScoresCorMatrix, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, bmassOutput$LogFile)[c("MergedDataSources", "MarginalSNPs", "ZScoresCorMatrix", "LogFile")]

	if (PrintMergedData == FALSE) {
		bmassOutput$MergedDataSources <- NULL
	}

	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Calculating logBFs and posteriors (if applicable).", sep=""), stderr())
	}
	
	bmassOutput[c("MarginalSNPs", "Models", "ModelPriors", "LogFile")] <- GetLogBFsFromData(DataSources, bmassOutput$MarginalSNPs, bmassOutput$ZScoresCorMatrix, SigmaAlphas, bmassOutput$LogFile)[c("MarginalSNPs", "Models", "ModelPriors", "LogFile")] 

	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "ModelPriors", "GWASlogBFMinThreshold", "LogFile")] <- DetermineAndApplyPriors(DataSources, bmassOutput$MarginalSNPs, GWASsnps, SigmaAlphas, bmassOutput$Models, bmassOutput$ModelPriors, ProvidedPriors, UseFlatPriors, GWASThreshFlag, GWASThreshValue, bmassSeedValue, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "ModelPriors", "GWASlogBFMinThreshold", "LogFile")]
	
	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Getting final list of MarginalSNPs, PreviousSNPs, and NewSNPs (where applicable).", sep=""), stderr())
	}

	bmassOutput[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")] <- FinalizeAndFormatResults(DataSources, bmassOutput$MarginalSNPs, bmassOutput$PreviousSNPs, GWASsnps, bmassOutput$GWASlogBFMinThreshold, SigmaAlphas, bmassOutput$Models, bmassOutput$ModelPriors, NminThreshold, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, bmassOutput$LogFile)[c("MarginalSNPs", "PreviousSNPs", "NewSNPs", "LogFile")]

	if (PrintProgress == TRUE) {
        	write(paste(format(Sys.time()), " -- Finishing bmass, exiting.", sep=""), stderr())
	}

	return(bmassOutput)
}
