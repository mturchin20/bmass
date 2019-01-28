#' Bayesian multivariate analysis of summary statistics (bmass)
#'
#' Run bmass on a set of phenotypes that each have univariate GWAS statistics on the same set of SNPs 
#' 
#' @usage bmass <- function (DataSources, GWASsnps = NULL, GWASThreshFlag = FALSE, GWASThreshValue = 5e-8, NminThreshold = 0)
#"
#' @param DataSources A string indicating the variable names of the input datafiles and phenotypes 
#' @param GWASsnps A data.table containing rows of SNPs that were univariate genome-wide significant in the phenotypes being used for analysis; GWASsnps input file should have two columns, one for chromosome and another for basepair position (with column headers of "Chr" and "BP")
#' @param GWASThreshFlag A logical TRUE/FALSE flag that indicates whether to threshold input GWASsnps list by a univariate GWAS p-value or not (eg the input GWASsnps list contains variants that are significant from discovery + replication data, but the input summary statistics are just from the discovery cohort). Default is TRUE. 
#' @param GWASThreshValue A numerical value indicating the univariate p-value threshold to use in conjunction with the GWASThreshFlag. Default is 5e-8. 
#' @param SNPMarginalUnivariateThreshold A numerical value indicating the univariate p-value threshold to use when collecting marginally significant SNPs for final bmass analysis. Default is 1e-6.  
#' @param SNPMarginalUnivariateThreshold A numerical value indicating the basic multivariate p-value threshold to use when collecting marginally significant SNPs for final bmasss analysis. Default is 1e-6.
#' @param NminThreshold A numerical value that indicates a sample size threshold to use where SNPs below which are removed. Default is 0.
#' @param PrintMergedData A logical TRUE/FALSE flag that indicates whether the intermediary 'merged datafile' should be included in the final bmass output; this file combines all the phenotypes for every SNP provided just prior to thresholding for marginally significant SNPs. Default is FALSE. 
#' @param PrintProgress A logical TRUE/FALSE flag that indicates whether progress statements should be printed to stderr() during the course of running bmass() or not. Default is FALSE.
#'
#' @return A list containing model, SNP, and posterior information for both the previously significant univariate SNPs ("PreviousSNPs") and the newly significant multivariate SNPs ("NewSNPs"). For a full breakdown of the bmass() output list structure, please see the associated README and/or vignettes. 
#'
#' @examples
#' bmass(c("HDL", "LDL", "TG", "TC"), GWASsnps, NminThreshold = 50000) 
#' bmass(c("HDL", "LDL", "TG", "TC"), GWASsnps, GWASThreshValue = 1e-8, NminThreshold = 50000, PrintProgress=TRUE) 
#' bmass(c("HDL", "LDL", "TG", "TC"), GWASsnps, GWASThreshFlag = FALSE, SNPMarginalUnivariateThreshold = 1e-4, SNPMarginalMultivariateThreshold = 1e-4, PrintMergedData = TRUE) 
#' bmassOutput <- list(); bmassOutput <- bmass(c("HDL", "LDL", "TG", "TC"), GWASsnps, NminThreshold = 50000) 
#'
bmass <- function (DataSources, GWASsnps = NULL, MergedDataSources = NULL, ZScoresCorMatrix = NULL, ExpectedColumnNames = c("Chr", "BP", "Marker", "MAF", "A1", "Direction", "pValue", "N"), GWASsnps_AnnotateWindow = 5e5, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), ProvidedPriors=NULL, UseFlatPriors=FALSE, GWASThreshFlag = TRUE, GWASThreshValue = 5e-8, NminThreshold = 0, PruneMarginalSNPs=TRUE, PruneMarginalSNPs_bpWindow=5e5, PrintMergedData=FALSE, PrintProgress=FALSE, bmassSeedValue=1) {

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
#		bmassOutput$LogFile <- CheckMergedDataSources(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, bmassOutput$MergedDataSources, ProvidedPriors, UseFlatPriors, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, bmassOutput$LogFile)
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
        	write(paste(format(Sys.time()), " -- finishing bmass, exiting.", sep=""), stderr())
	}

	return(bmassOutput)

}

