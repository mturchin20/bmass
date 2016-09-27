#' Checking and preparing input datafiles.
#' 
#' Description 
#' 
#' @param x Something.
#' @param y Something2.
####' @param DataFileList Comma-separated list of all the datafile names being analyzed (these are assumed to be the names of the phenotypes being analyzed). This should be in the same order as DataFileLocations. The default value for this paramter is NULL. If NULL, the list of datafile names are derived from the list of datafile locations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList. The default value for this parameter is NULL. Users must supply a list of datafile locations to run bmass.
#' @param ExpectedColumnNames Comma-separated list of the expected column names to be found at the head of each datafile. The default value for this parameter is Chr,BP,A1,MAF,Direction,pValue,N. Users should not supply or alter this parameter.
#' @param OutputFileBase
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)

#Data1 <- read.table("../data/TestData1.txt", header=T)
#Data2 <- read.table("../data/TestData2.txt", header=T)
#SigSNPs <- read.table("../data/TestData1.GWASsnps.txt", header=T)
##source("PrepareData.R")
#library("devtools")
#devtools::load_all()
#bmass(c("Data1", "Data2"), GWASsnps=SigSNPs, NminThreshold = 2000, bmassSeedValue=NULL)
#bmassOutput1 <- bmass(c("Data1", "Data2"), GWASsnps=SigSNPs, bmassSeedValue=NULL)
#ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N")
#DataList <- c("Data1", "Data2")

CheckCharacterFormat <- function (DataSource1) {
	returnValue <- FALSE
	
	if (is.character(DataSource1)) {
		returnValue <- TRUE
	}

	return(returnValue)
}

CheckVariableExists <- function (DataSource1) {
	returnValue <- FALSE
	
	if (exists(DataSource1)) {
		returnValue <- TRUE
	}

	return(returnValue)
}

CheckDataFrameFormat <- function (DataSource1) {
	returnValue <- FALSE
	
	if (is.data.frame(eval(parse(text=DataSource1)))) {
			returnValue <- TRUE
	}

	return(returnValue)
}

CheckDataSourceHeaders <- function (DataSources1, ExpectedColumnNames1) {
	returnValueVector <- c()
	for (DataSource1 in DataSources1) {
		returnValue <- TRUE
		DataSource1ColumnNames <- names(eval(parse(text=DataSource1)))
		for (ExpectedColumnName in ExpectedColumnNames1) {
			if (! ExpectedColumnName %in% DataSource1ColumnNames) {
				returnValue <- FALSE
			}
		}
		returnValueVector <- c(returnValueVector, returnValue)
	}

	return(returnValueVector)
}

CheckDataSourceDirectionColumn <- function (DataSources1) {
	returnValue <- TRUE
	for (DirectionValue in eval(parse(text=paste(DataSources1, "$Direction", sep="")))) {
		if ((DirectionValue != "+") && (DirectionValue != "-")) {
			returnValue <- FALSE
		}
	}
	return(returnValue)
}

#~~~
#> txt3 <- "nana,nana2"
#> grep(",", txt3)
#[1] 1
#> ?grep
#> txt3 <- "nana,nana2,nana3"
#> grep(",", txt3)
#[1] 1
#> ?grep
#> grep(",", c(txt3,txt3))
#[1] 1 2
#> grep(",", c(txt3,txt3,txt3))
#[1] 1 2 3
#> grep(",", c(txt3,txt3,txt3,txt2))
#[1] 1 2 3
#> grep(",", c(txt3,txt3,txt2,txt3))
#[1] 1 2 4
#~~~
#~~~
#> txt1 <- "nana"
#> txt2 <- "nana/nan2"
#> strsplit(txt1, "/")
#[[1]]
#[1] "nana"
#
#> strsplit(txt2, "/")
#[[1]]
#[1] "nana" "nan2"
#
#> strsplit(txt2, "/")[[1]][1]
#[1] "nana"
#> strsplit(txt2, "/")[[1]][2]
#[1] "nan2"
#> length(strsplit(txt2, "/")[[1]])  
#[1] 2
#> strsplit(txt2, "/")[[1]][length(strsplit(txt2, "/")[[1]])]
#[1] "nan2"
#~~~

#This is going to be the main function that goes through each of the steps from beginning to end. Hypothetically, all the other functions presented here should be used through the PrepareData process (or as a subfunction of one of the functions being used in PrepareData)
#PrepareData <- function (ExpectedColumnNames, DataFileLocations, OutputFileBase) { #20160814 NOTE -- Changing direction and just assuming input is a single vector that contains all the proper data.frame datasources and working from there. Final output will be a list that has all the output. 'Logfile' will just be a variable included in final list output, developed by continula 'rbind' calls with text output additions. Also deciding to move 'PrepareData' to just a 'MainWorkFlow' or 'Main' that I'll dev in each sub R package and then eventually move to a main source.  
#bmass <- function (DataSources, GWASsnps=NULL, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), MergedDataSources=NULL, ProvidedPriors=NULL, UseFlatPriors=FALSE, PruneMarginalHits=TRUE, PruneMarginalHits_bpWindow=5e5, SNPMarginalUnivariateThreshold = 1e-6, SNPMarginalMultivariateThreshold = 1e-6, NminThreshold = 0, bmassSeedValue=NULL) {
CheckIndividualDataSources <- function (DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile) {

	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- beginning bmass.", sep=""))

	#20160823 CHECK_0: Prob -- list of Matthew functions specifically to double-check, go through, go over
	#	collapse
	#	em.priorprobs

	#Loading and checking data
	#~~~~~~

	if (is.null(MergedDataSources)) {

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- beginning DataSources checks.", sep=""))
	
		####Candidate For Unit Tests####
		if (!is.vector(DataSources)) {
			stop(Sys.time(), " -- input variable DataSources not in vector format. bmass expects DataSources to be a vector of strings. Please fix and rerun bmass.") 
		}
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed vector check.", sep=""))

		####Candidate For Unit Tests####
		DataSourcesCheckCharacters <- sapply(DataSources, CheckCharacterFormat)	
		if (FALSE %in% DataSourcesCheckCharacters) {
			stop(Sys.time(), " -- the following entries in DataSources were not found as characters. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckCharacters])
		}
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed string check.", sep=""))

		####Candidate For Unit Tests####
		DataSourcesCheckExists <- sapply(DataSources, CheckVariableExists)
		if (FALSE %in% DataSourcesCheckExists) {
			stop(Sys.time(), " -- the variables associated with the following entries in DataSources were not found to exist. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckExists])
		}

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed exists check.", sep=""))
		
		#20160814 20160814 CHECK_1 -- Prob: Create way to specify which files are causing the data.frame failure here? Maybe change DataFrameCheckValues into a vector of true/false statements and then convert the trues to their text names as posible outputs? Soln: Moved to a format where returning TRUE and FALSE statements in a vector, and then pass that vector to DataSources character vector to get proper output. 
		####Candidate For Unit Tests####
		DataSourcesCheckDataFrames <- sapply(DataSources, CheckDataFrameFormat)
		if (FALSE %in% DataSourcesCheckDataFrames) {
			stop(Sys.time(), " -- the following data sources are not formatted as data.frames. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckDataFrames])
		}

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed data.frame check.", sep=""))

		####Candidate For Unit Tests####
		DataSourcesCheckHeaderNames <- CheckDataSourceHeaders(DataSources, ExpectedColumnNames)
		if (FALSE %in% DataSourcesCheckHeaderNames) {
			stop(Sys.time(), " -- the following data sources do not have all the expected column headers. The expected column headers are \"", paste(ExpectedColumnNames, collapse=" "), "\". Please fix and rerun bmass: ", DataSources[!DataSourcesCheckHeaderNames])
		}

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed column headers check.", sep=""))
		
		
		#20160901 CHECK_0 -- Prob: Do X/23 chr column conversion stuff first here
		###### do thissssss
		#
		#Check for 'X' chromosomes and convert to 23 & print warning/logfile statement
		#
		###### do thissssss
	
		#20160901 CHECK_0 -- Prob: Check multiple columns in entries to make sure input types/classes are as what's exepected? Eg chr/bp numeric, alleles chars, N and pvals numeric? Do X/23 conversion stuff before ondividual data column class/type checks first
		####Candidate For Unit Tests####
		DataSourcesCheckDirectionColumn <- sapply(DataSources, CheckDataSourceDirectionColumn)
		if (FALSE %in% DataSourcesCheckDirectionColumn) {
			stop(Sys.time(), " -- the following data sources have entries other than + and - in the Direction column. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckDirectionColumn])
		}
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed Direction column check.", sep=""))

		###### do thissssss below
		#
		#Check that if ProvidedPriors is non-NULL then it has the correct length and that all entries are numeric. 
		#
		###### do thissssss below
		####Candidate For Unit Tests####
		if (!is.null(ProvidedPriors)) {
			if (!is.vector(ProvidedPriors)) {
				stop(Sys.time(), " -- ProvidedPriors input is not in vector format. Please fix and rerun bmass.")
			}
			#20160901 CHECK_0 -- Prob: How stringent should test types be for input variables? Eg need to be testing ProvidedPriors is numeric, and other specific input variable classes too?
			if (!is.numeric(ProvidedPriors)) {
				stop(Sys.time(), " -- ProvidedPriors input is returning false for is.numeric(). Please ensure all entries are numeric and then rerun bmass.")
			}
			if (length(ProvidedPriors) != 3^length(DataSources)) { 
				stop(Sys.time(), " -- The number of entries in ProvidedPriors does not equal 3 ^ the number of datasets passed to DataSources (ie 3 ^ ", as.character(length(DataSources)), " = ", as.character(3^length(DataSources)), "). Please fix and rerun bmass.")
			}
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors was provided and passed checks.", sep=""))
		}
	
		#20160902 CHECK_0 -- Prob: Go through all other input variables and make sure they are the expected formats/inputs, eg numeric, character, etcetc
		#bmass <- function (DataSources, GWASsnps=NULL, SNPMarginalpValThreshold=1e-6, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), MergedDataSources=NULL, ProvidedPriors=NULL, UseFlatPriors=FALSE, PruneMarginalHits=TRUE, PruneMarginalHits_bpWindow=5e5, NminThreshold = 0, bmassSeedValue=NULL) 

		if (!is.numeric(NminThreshold)) {
			stop(Sys.time(), " -- .")
		}

	}
	else {
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- MergedDataSources file provided, going through data checks.", sep=""))

		#20160901 CHECK_0 -- Prob: Get this done at some point soon'ish, not soon after first rough draft completed for everything
		###### do thissssss
		#
		#Routine for checks if MergedDataSources file provided?
		#
		###### do thissssss

	}













	if (!is.null(ProvidedPriors)) {
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
		MarginalHits_logBFs$prior <- ProvidedPriors
	}

#	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up priors.", sep=""))

#	GWASHits_EBprior <- NULL
#	FlatUnif_EBprior <- NULL
#	logBF_min	
	MarginalHits_logBFs_Stacked_AvgwPrior <- NULL
	MarginalHits_logBFs_Stacked_AvgwPrior_Min <- NULL
	Priors_Used <- NULL
        if (!is.null(ProvidedPriors)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
                MarginalHits_logBFs$prior <- ProvidedPriors
        	Priors_Used <- ProvidedPriors
	} else if (is.null(GWASsnps) || UseFlatPriors == TRUE) {
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up flat-tiered priors, GWASnps either not provided or flat prior explicitly requested.", sep=""))
	
		Prior_FlatUnif <- normalize(rep(c(0,MarginalHits_logBFs$prior[-1]),length(SigmaAlphas)))
		#nsigma=length(sigmaa)
		#origprior = rep(c(0,lbf$prior[-1]),nsigma)
		#origprior = normalize(origprior)
		
		#Prior_FlatUnif 
	
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and flat-tiered priors.", sep=""))
	
		MarginalHits_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalHits_logBFs_Stacked, Prior_FlatUnif) 
		Priors_Used <- Prior_FlatUnif	
	
		#lbf.av.origprior.glhits = lbf.av(lbf.glhits,origprior)

		#Add summary stats to marginal SNPs
		#Add SNPs x Model matrix with prior*logBFs as entries, summed (or avg'd??) across SigmaAlphas
	
	}
	else {
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up GWAS trained priors and analyzing GWAS hits since GWASsnps provided.", sep=""))
		PreviousSNPs <- MarginalHits[MarginalHits$GWASannot==1,]
		PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalHits_logBFs_Stacked[,MarginalHits$GWASannot==1]) #Matrix of nSigmaAlphas x nSNPs
		
		#20160822 20160823 CHECK_1 -- Prob: Do GWAS hit analysis/work here Soln: Wrote up a few first-level GWAS hit results to get things started. Certainly will undergo further revisions down the line but fine strating point for now.
		
		Prior_PreviousSNPsEB <- em.priorprobs(PreviousSNPs_logBFs_Stacked, MarginalHits_logBFs$prior, 100) #Vector with nModels*nSigmaAlphas entries
		#20160823 CHECK_0: Prob -- double check use of em.priorprobs here too with runif prior starting point
		#20160823 CHECK_0: Prob -- Do multiple runs of em.priorprobs and figure out way to compare them for consistency?
		Prior_PreviousSNPsEB_check2 <- em.priorprobs(PreviousSNPs_logBFs_Stacked, MarginalHits_logBFs$prior*runif(length(MarginalHits_logBFs$prior)), 100)
	
		###### do thissssss
		#
		# Run multiple EMs to check/test for convergence? 
		#
		###### do thissssss

		Prior_PreviousSNPsEB_Collapsed <- CollapseSigmaAlphasTogether(Prior_PreviousSNPsEB, length(SigmaAlphas))
		Prior_PreviousSNPsEB_check2_Collapsed <- CollapseSigmaAlphasTogether(Prior_PreviousSNPsEB_check2, length(SigmaAlphas))

		PreviousSNPs_PosteriorProbs <- posteriorprob(PreviousSNPs_logBFs_Stacked, Prior_PreviousSNPsEB) #Matrix of nModels*nSigmaAlphas x nSNPs 
		PreviousSNPs_PosteriorProbs_Collapsed <- apply(PreviousSNPs_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
#		print(PreviousSNPs_PosteriorProbs)
#		print(PreviousSNPs_PosteriorProbs_Collapsed)
		
###		PrevHits$BestModels <- apply(MarginalHits_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")
		PreviousSNPs$BestModel <- apply(MarginalHits_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")	
		PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, max)
#		print(PrevHits$BestModels)
		print(PreviousSNPs)

		## this returns a list with elements  pU pD and pI
		#marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

		Prior_PreviousSNPsEB_ModelMatrix <- cbind(MarginalHits_logBFs$gamma, Prior_PreviousSNPsEB_Collapsed)[order(Prior_PreviousSNPsEB_Collapsed, decreasing=TRUE),]
		Prior_PreviousSNPsEB_ModelMatrix <- data.frame(cbind(Prior_PreviousSNPsEB_ModelMatrix, cumsum(Prior_PreviousSNPsEB_ModelMatrix[,ncol(Prior_PreviousSNPsEB_ModelMatrix)])))
		colnames(Prior_PreviousSNPsEB_ModelMatrix) <- c(DataSources, "pValue", "Cumm_pValue")
		print(Prior_PreviousSNPsEB_ModelMatrix)

		Prior_PreviousSNPsEB_ModelMatrix_AllAssoc <- apply((Prior_PreviousSNPsEB_ModelMatrix[,1:length(DataSources)]>0), 1, sum) == length(DataSources)
		Prior_PreviousSNPsEB_ModelMatrix_AllAssoc_pValSum <- sum(Prior_PreviousSNPsEB_ModelMatrix[Prior_PreviousSNPsEB_ModelMatrix_AllAssoc,ncol(Prior_PreviousSNPsEB_ModelMatrix)-1])
		Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc <- apply((Prior_PreviousSNPsEB_ModelMatrix[,1:length(DataSources)]>0), 1, sum) == length(DataSources)-1
		Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc_pValSum <- sum(Prior_PreviousSNPsEB_ModelMatrix[Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc,ncol(Prior_PreviousSNPsEB_ModelMatrix)-1])

		Prior_PreviousSNPsEB_ModelMatrix_pValSupport <- c(Prior_PreviousSNPsEB_ModelMatrix_AllAssoc_pValSum, Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc_pValSum)
		Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names <- c("AllAssoc", "AllBut1Assoc")
		for (DataSource in DataSources) {
			eval(parse(text=paste("Prior_PreviousSNPsEB_ModelMatrix_pValSupport <- c(Prior_PreviousSNPsEB_ModelMatrix_pValSupport, sum(Prior_PreviousSNPsEB_ModelMatrix[Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc & Prior_PreviousSNPsEB_ModelMatrix[,\"", DataSource, "\"] == 0, ncol(Prior_PreviousSNPsEB_ModelMatrix)-1]))", sep="")))
			Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names <- c(Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names, paste("AllBut", DataSource, "Assoc", sep=""))
		}
		names(Prior_PreviousSNPsEB_ModelMatrix_pValSupport) <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names

		#20160905 NOTE -- code below eventually moved to a downstream follow-up analysis helper function. Being kept here for the moment just while focus on rewriting and cleaning up more pertinent parts of the code first
		PrevHits <- list()
		PrevHits$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
#		PreviousSNPs$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
#		print(Prior_PreviousSNPsEB_ModelMatrix_pValSupport)
		print(PrevHits$ModelCategories_CummpValues)
	
###?		MarginalHits_logBFs_Stacked_AvgwPrior <- ()
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and trained priors.", sep=""))
		
		MarginalHits_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalHits_logBFs_Stacked, Prior_PreviousSNPsEB) 
		MarginalHits_logBFs_Stacked_AvgwPrior_Min <- min(MarginalHits_logBFs_Stacked_AvgwPrior)
		Priors_Used <- Prior_PreviousSNPsEB
	}

	#20160901 CHECK_0 -- Prob: Do something more substantive here? A better error message, or give just a warning instead? Don't exist program?
	####Candidate For Unit Tests####
	if (is.null(MarginalHits_logBFs_Stacked_AvgwPrior)) {
		stop(Sys.time(), " -- No average log BFs were returned from method. Check if all input variables are as the method expects.") 
	}

	MarginalHits$LogBFWeightedAvg <- MarginalHits_logBFs_Stacked_AvgwPrior



















	#Pruning marginal hits by LogBFWeightedAvg if requested
	if (PruneMarginalHits == TRUE) {
		#20160901 CHECK_0 -- Prob: Go over indepthits function, rewrite, or just lightly edit? redo names, double-check functionality? def get some unit testing in there
		#MarginalHits <- MarginalHits[indephits(MarginalHits$LogBFWeightedAvg, MarginalHits$Chr, MarginalHits$BP, T=PruneMarginalHits_bpWindow)==1,]
		MarginalHits_PrunedList <- indephits(MarginalHits$LogBFWeightedAvg, MarginalHits$Chr, MarginalHits$BP, T=PruneMarginalHits_bpWindow)
		MarginalHits <- MarginalHits[MarginalHits_PrunedList==1,]	
		MarginalHits_logBFs$lbf <- lapply(MarginalHits_logBFs$lbf, function(x) { return(x[,MarginalHits_PrunedList==1]) })
		MarginalHits_logBFs_Stacked <- MarginalHits_logBFs_Stacked[,MarginalHits_PrunedList==1]	
	}

	print(MarginalHits)
#	print(summary(MarginalHits_logBFs$lbf))

	#Preparing log Bayes factors matrix, Gammas x SNPs	
	#bmassOutput$logBFs <- GetSumAcrossSigmaAlphas_withPriors(MarginalHits_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalHits_logBFs_Stacked)), ncol=ncol(MarginalHits_logBFs_Stacked), byrow=FALSE), nrow(MarginalHits_logBFs$gamma), length(SigmaAlphas))
	MarginalHits_logBFs_Stacked_SigmaAlphasSummed <- GetSumAcrossSigmaAlphas_withPriors(MarginalHits_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalHits_logBFs_Stacked)), ncol=ncol(MarginalHits_logBFs_Stacked), byrow=FALSE), nrow(MarginalHits_logBFs$gamma), length(SigmaAlphas))
	MarginalHits_logBFs_Stacked_SigmaAlphasSummed <- cbind(MarginalHits_logBFs$gamma, MarginalHits_logBFs_Stacked_SigmaAlphasSummed)
	colnames(MarginalHits_logBFs_Stacked_SigmaAlphasSummed) <- c(DataSources, MarginalHits$ChrBP)

	print(MarginalHits_logBFs_Stacked_SigmaAlphasSummed)
	print(Priors_Used)
	print(MarginalHits_logBFs_Stacked)

	#Preparing posterior probabilities, Gammas x SNPs
	MarginalHits_logBFs_Stacked_PosteriorProbabilities <- posteriorprob(MarginalHits_logBFs_Stacked, Priors_Used)
	MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed <- apply(MarginalHits_logBFs_Stacked_PosteriorProbabilities, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)
)
	MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed <- cbind(MarginalHits_logBFs$gamma, MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed)
	colnames(MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed) <- c(DataSources, MarginalHits$ChrBP)
	#MarginalHits_PosteriorProbabilities <- apply(posteriorprob(MarginalHits_logBFs_Stacked, Priors_Used), 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas))
	
	print(MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed)
	
	#20160905 CHECK_0 -- Prob: Convert either 'MarginalHits_logBFs_Stacked_PosteriorProbabilities' to 'PosteriorProbs' or 'PreviousSNPs_PosteriorProbs' to 'PreviousSNPs_PosteriorProbabilities'
	#20160905 CHECK_0 -- Prob: Change 'MarginalHits' to 'MarginalSNPs' probably?
	#20160905 NOTE -- Below code was not tested yet, just typing it out now to be included/tested later. This was something I wanted to start including anyways. 
	#MarginalHits$BestModel <- apply(MarginalHits_logBFs$gamma[apply(MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed, 2, which.max),], 1, paste, collapse="_")
	#MarginalHits$BestModel_Posterior <- apply(MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed, 2, max)
	
	#PreviousSNPs_PosteriorProbs <- posteriorprob(PreviousSNPs_logBFs_Stacked, Prior_PreviousSNPsEB) #Matrix of nModels*nSigmaAlphas x nSNPs 
	#PreviousSNPs_PosteriorProbs_Collapsed <- apply(PreviousSNPs_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
	#PreviousSNPs$BestModel <- apply(MarginalHits_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")	
	#PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, max)
	
	#lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 81, 14)
	#lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
	#lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
	#lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]

	#20160905 CHECK_0 -- Prob: Move this below intialization section to proper beginning of .R file code as necessary once change/reorganization occurs
	NewSNPs <- NULL	

	print(MarginalHits_logBFs_Stacked_AvgwPrior_Min)
	
	#Determining new hits if GWASsnps were provided to determine minimum MarginalHits_logBFs_Stacked_AvgwPrior value threshold	
	#20160902 CHECK_0 -- Prob: Check that MarginalHits_logBFs_Stacked_AvgwPrior_Min is non null here?
	if (!is.null(GWASsnps)) {
		if (is.null(MarginalHits_logBFs_Stacked_AvgwPrior_Min)) {
			stop(Sys.time(), " -- MarginalHits_logBFs_Stacked_AvgwPrior_Min is NULL despite GWASsnps being provided. Unexpected error.")
		}
		NewSNPs <- MarginalHits[MarginalHits$GWASannot == 0 & MarginalHits$LogBFWeightedAvg >= MarginalHits_logBFs_Stacked_AvgwPrior_Min & MarginalHits$Nmin >= NminThreshold,]
	}
	

	#Preparing final return variable bmassOutput
	bmassOutput$MarginalSNPs$SNPs <- MarginalHits
	print(dim(MarginalHits_logBFs_Stacked))
	print(length(Priors_Used))
	bmassOutput$ModelPriors <- Priors_Used
	####Candidate For Unit Tests####
	print(dim(matrix(rep(Priors_Used, ncol(MarginalHits_logBFs_Stacked)), ncol=ncol(MarginalHits_logBFs_Stacked), byrow=FALSE)))
	#bmassOutput$logBFs <- GetSumAcrossSigmaAlphas_withPriors(MarginalHits_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalHits_logBFs_Stacked)), ncol=ncol(MarginalHits_logBFs_Stacked), byrow=FALSE), nrow(MarginalHits_logBFs$gamma), length(SigmaAlphas)) 
	bmassOutput$MarginalSNPs$logBFs <- MarginalHits_logBFs_Stacked_SigmaAlphasSummed 
	print(dim(bmassOutput$logBFs))
	print(bmassOutput$logBFs)
	bmassOutput$MarginalSNPs$Posteriors <- MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed
	bmassOutput$GWASlogBFMinThreshold <- MarginalHits_logBFs_Stacked_AvgwPrior_Min
	bmassOutput$NewSNPs$SNPs <- NewSNPs
	print(MarginalHits_logBFs_Stacked_AvgwPrior_Min)
	print(bmassOutput$NewSNPs)
		
	#PreviousSNPs$BestModel <- apply(MarginalHits_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")	
	#PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, max)
	
	#bmassOutput <- list()
	#bmassOutput$ModelPriors <- NULL
	#bmassOutput$MarginalSNPs <- list() 
	#bmassOutput$MarginalSNPs$logBFs <- NULL
	#bmassOutput$MarginalSNPs$Posteriors <- NULL
	#bmassOutput$NewSNPs <- list()
	#bmassOutput$NewSNPs$SNPs <- NULL
	#bmassOutput$NewSNPs$logBFs <- NULL
	#bmassOutput$NewSNPs$Posteriors <- NULL
	#bmassOutput$GWASlogBFMinThreshold <- NULL
	#bmassOutput$PreviousSNPs <- list()
	#bmassOutput$PreviousSNPs$SNPs <- NULL
	#bmassOutput$PreviousSNPs$logBFs <- NULL
	#bmassOutput$PreviousSNPs$Posteriors <- NULL


	#Returning final output variable, bmassOutput

	#return(bmassOutput)

	##extract lbfs for all the new hits
	#lbf.newhits= lbf.bigmat[,gl$nmin>20000]
	#lbf.newhits= lbf.newhits[,l==1]
	#lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>5.083439 & sub$nmin>20000]
	#pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
	#pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))
	#
	#ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
	#pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
	#bestclass= apply(pp.newhits.classmatrix, 2,which.max)
	#cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

}











#Expectation
#Vector input of data.frame data sources 

#' @param DataSources input names of datasources
#' @param LogFile rbind'ing vector text output
#' @param ExpectedColumnNames Comma-separ

#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ

#if (FALSE) {
#~~~
#TimeTag 20160902 -- Right above is dealing with looking at pruning subsection output being correct, below is having to do with NewHits section and making sure pruning is extending into the expected places of logBFs and PosteriorProbabilities 
#      Length Class  Mode   
# [1,] 45     -none- numeric
# [2,] 45     -none- numeric
# [3,] 45     -none- numeric
# [4,] 45     -none- numeric
# [5,] 45     -none- numeric
# [6,] 45     -none- numeric
# [7,] 45     -none- numeric
# [8,] 45     -none- numeric
# [9,] 45     -none- numeric
#[10,] 45     -none- numeric
#[11,] 45     -none- numeric
#[12,] 45     -none- numeric
#[13,] 45     -none- numeric
#[14,] 45     -none- numeric
#      Length Class  Mode   
# [1,] 36     -none- numeric
# [2,] 36     -none- numeric
# [3,] 36     -none- numeric
# [4,] 36     -none- numeric
# [5,] 36     -none- numeric
# [6,] 36     -none- numeric
# [7,] 36     -none- numeric
# [8,] 36     -none- numeric
# [9,] 36     -none- numeric
#[10,] 36     -none- numeric
#[11,] 36     -none- numeric
#[12,] 36     -none- numeric
#[13,] 36     -none- numeric
#[14,] 36     -none- numeric
#[1] 126   5
#[1] 126
#[1] 126   5
#[1] 9 5
#             [,1]        [,2]        [,3]         [,4]       [,5]
# [1,]    0.000000    0.000000    0.000000    0.0000000    0.00000
# [2,]    1.952878    1.331036    1.130067    1.0750168   34.86177
# [3,]    0.000000    0.000000    0.000000    0.0000000    0.00000
# [4,]    2.405875    1.139964    0.626045    0.8584684   38.74820
# [5,]  -54.397661  -59.826155    5.711266  -60.1456968  -24.36255
# [6,] -213.542751 -218.251940 -218.979761 -218.3956795 -213.78594
# [7,]    0.000000    0.000000    0.000000    0.0000000    0.00000
# [8,] -207.433088 -211.552755 -211.904431 -211.6312350 -210.84231
# [9,]    0.000000    0.000000    0.000000    0.0000000    0.00000
#[1] 126   4
#[1] 126
#[1] 126   4
#[1] 9 4
#             [,1]        [,2]         [,3]       [,4]
# [1,]    0.000000    0.000000    0.0000000    0.00000
# [2,]    1.952878    1.331036    1.0750168   34.86177
# [3,]    0.000000    0.000000    0.0000000    0.00000
# [4,]    2.405875    1.139964    0.8584684   38.74820
# [5,]  -54.397661  -59.826155  -60.1456968  -24.36255
# [6,] -213.542751 -218.251940 -218.3956795 -213.78594
# [7,]    0.000000    0.000000    0.0000000    0.00000
# [8,] -207.433088 -211.552755 -211.6312350 -210.84231
# [9,]    0.000000    0.000000    0.0000000    0.00000
#[1] 5.711266
#    ChrBP Chr    BP A1  MAF Data1_Direction Data1_pValue Data1_N Data1_ZScore
#6 3_21000   3 21000  G 0.37               +        5e-08    2582     5.451310
#9  4_7000   4  7000  G 0.15               +        0e+00    2514     5.931598
#  Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot    mvstat
#6               +        9e-08    2510     5.345837         0  36.50763
#9               -        0e+00    2514    -7.348796         0 219.57617
#  mvstat_log10pVal  unistat unistat_log10pVal Nmin LogBFWeightedAvg
#6         7.927532 29.71679           7.30103 2510         6.257912
#9        47.680359 54.00480          12.69897 2514        45.346014
#}





