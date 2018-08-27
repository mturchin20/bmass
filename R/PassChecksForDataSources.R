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

CheckCharacterClass <- function (DataSource1) {
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

#20160929 CHECK_0: Prob -- Come back and check new functionality version below when running on 2nd version of bmass
CheckDataSourceDirectionColumn <- function (DataSources1) {
	returnValue <- TRUE
	if (length(eval(parse(text=paste(DataSources1, "$Direction", sep="")))[eval(parse(text=paste(DataSources1, "$Direction", sep=""))) != "+" & eval(parse(text=paste(DataSources1, "$Direction", sep=""))) != "-"]) > 0) {
		returnValue <- FALSE
	}
	return(returnValue)
}

CheckDataSourceMAFIsMAF <- function (DataSources1) {
	returnValue <- TRUE
	if (TRUE %in% eval(parse(text=paste(DataSources1, "$MAF > .5 ", sep="")))) {
		returnValue <- FALSE
	}
	return(returnValue)
}

CheckDataSourceMAFFixed <- function (DataSources1) {
	returnValue <- TRUE
	returnVector <- eval(parse(text=paste(DataSources1, "$MAF == 0 | ", DataSources1, "$MAF == 1 ", sep="")))
#	write(head(returnVector), stderr());
	if (TRUE %in% returnVector) {
		returnValue <- FALSE
#		write(head(eval(parse(text=paste(DataSources1, " <- ", DataSources1, "[", !returnVector, ",]", sep="")))), stderr())
#		eval(parse(text=paste(DataSources1, " <- ", DataSources1, "[", returnVector, ",]", sep=""))) 
	}
	return(returnValue)
}

CheckIndividualDataSources <- function (DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile) {

	#20160823 CHECK_0: Prob -- list of Matthew functions specifically to double-check, go through, go over
	#	collapse
	#	em.priorprobs

	#Loading and checking data
	#~~~~~~

	if (is.null(MergedDataSources)) {

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- beginning DataSources checks.", sep=""))
	
		if (!is.vector(DataSources)) {
			stop(Sys.time(), " -- input variable DataSources not in vector format. bmass expects DataSources to be a vector of strings. Please fix and rerun bmass.") 
		}
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed vector check.", sep=""))

		DataSourcesCheckCharacters <- sapply(DataSources, CheckCharacterClass)	
		if (FALSE %in% DataSourcesCheckCharacters) {
			stop(Sys.time(), " -- the following entries in DataSources were not found as characters. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckCharacters], collapse=" "))
		}
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed string check.", sep=""))

		DataSourcesCheckExists <- sapply(DataSources, CheckVariableExists)
		if (FALSE %in% DataSourcesCheckExists) {
			stop(Sys.time(), " -- the variables associated with the following entries in DataSources were not found to exist. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckExists], collapse=" "))
		}

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed exists check.", sep=""))
		
		#20160814 20160814 CHECK_1 -- Prob: Create way to specify which files are causing the data.frame failure here? Maybe change DataFrameCheckValues into a vector of true/false statements and then convert the trues to their text names as posible outputs? Soln: Moved to a format where returning TRUE and FALSE statements in a vector, and then pass that vector to DataSources character vector to get proper output. 
		DataSourcesCheckDataFrames <- sapply(DataSources, CheckDataFrameFormat)
		if (FALSE %in% DataSourcesCheckDataFrames) {
			stop(Sys.time(), " -- the following data sources are not formatted as data.frames. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckDataFrames], collapse=" "))
		}

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed data.frame check.", sep=""))

		DataSourcesCheckHeaderNames <- CheckDataSourceHeaders(DataSources, ExpectedColumnNames)
		if (FALSE %in% DataSourcesCheckHeaderNames) {
			stop(Sys.time(), " -- the following data sources do not have all the expected column headers. The expected column headers are \"", paste(ExpectedColumnNames, collapse=" "), "\". Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckHeaderNames], collapse=" "))
		}

		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed column headers check.", sep=""))
		
		#20160901 CHECK_0 -- Prob: Do X/23 chr column conversion stuff first here
		###### do thissssss
		#
		#Check for 'X' chromosomes and convert to 23 & print warning/logfile statement
		#
		###### do thissssss
	
		#20160901 CHECK_0 -- Prob: Check multiple columns in entries to make sure input types/classes are as what's exepected? Eg chr/bp numeric, alleles chars, N and pvals numeric? Do X/23 conversion stuff before ondividual data column class/type checks first
		DataSourcesCheckDirectionColumn <- sapply(DataSources, CheckDataSourceDirectionColumn)
		if (FALSE %in% DataSourcesCheckDirectionColumn) {
			stop(Sys.time(), " -- the following data sources have entries other than + and - in the Direction column. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckDirectionColumn], collapse=" "))
		}
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed Direction column check.", sep=""))

		#20171018 CHECK_0 -- Prob: Give first-pass check of the below function
		DataSourcesCheckMAFIsMAF <- sapply(DataSources, CheckDataSourceMAFIsMAF)
		if (FALSE %in% DataSourcesCheckMAFIsMAF) {
			stop(Sys.time(), " -- the following data sources have variants whose MAF entry are > .5; bmass expects the MAF column to only have values <= .5. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckMAFIsMAF], collapse=" "))
		}
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed MAF column check.", sep=""))

		#20171018 CHECK_0 -- Prob: Give first-pass check of the below function
		DataSourcesCheckMAFFixed <- sapply(DataSources, CheckDataSourceMAFFixed) 
		if (FALSE %in% DataSourcesCheckMAFFixed) {
			stop(Sys.time(), " -- the following data sources have variants whose MAF are == 0 (or == 1); bmass expects only segregating variants (eg not fixed). Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckMAFFixed], collapse=" "))
		} else {
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed MAF fixed check.", sep=""))
		}

		#20160901 CHECK_0 -- Prob: How stringent should test types be for input variables? Eg need to be testing ProvidedPriors is numeric, and other specific input variable classes too?
		#20160930 CHECK_0 -- Prob: Go over this below section and make sure it contains everything wanted for when 'ProvidedPriors' is in fact provided?
		if (!is.null(ProvidedPriors)) {
			if (!is.vector(ProvidedPriors)) {
				stop(Sys.time(), " -- ProvidedPriors input is not in vector format. Please fix and rerun bmass.")
			}
			if (!is.numeric(ProvidedPriors)) {
				stop(Sys.time(), " -- ProvidedPriors input is returning false for is.numeric(). Please ensure all entries are numeric and then rerun bmass.")
			}
			#20170218 20171017 CHECK_1 -- Prob: come back here and check if 14 was overall needed or just a stopgap while doing some bug-testing. Originally there was no `*14`, it was just `!= 3^length(DataSources)) {` Soln: 14 corresponds to `length(SigmaAlphas)` which makes sense since there are priors not only for all models but all models across all `SigmaAlphas`
			if (length(ProvidedPriors) != 3^(length(DataSources)*length(SigmaAlphas))) { # 20171017 NOTE -- originally just had `*14` in place of where `length(SigmaAlphas)` is now
				stop(Sys.time(), " -- The number of entries in ProvidedPriors does not equal 3 ^ (the number of datasets passed to DataSources * length(SigmaAlphas)) (ie 3 ^ (", as.character(length(DataSources)), " * ", as.character(length(SigmaAlphas)), ") = ", as.character(3^(length(DataSources)*length(SigmaAlphas))), "). Please fix and rerun bmass.")
			}
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors was provided and passed checks.", sep=""))
		}
	
		#20160902 CHECK_0 -- Prob: Go through all other input variables and make sure they are the expected formats/inputs, eg numeric, character, etcetc
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

	return(LogFile)
}



#Expectation
#Vector input of data.frame data sources 

#' @param DataSources input names of datasources
#' @param LogFile rbind'ing vector text output
#' @param ExpectedColumnNames Comma-separ

#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ



