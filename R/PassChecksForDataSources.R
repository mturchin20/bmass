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





