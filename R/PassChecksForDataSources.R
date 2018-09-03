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
	if (TRUE %in% returnVector) {
		returnValue <- FALSE
	}
	return(returnValue)
}

CheckIndividualDataSources <- function (DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile) {

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
		
		#20160901 20180830 CHECK_10 -- Prob: Check multiple columns in entries to make sure input types/classes are as what's exepected? Eg chr/bp numeric, alleles chars, N and pvals numeric? Do X/23 conversion stuff before ondividual data column class/type checks first Soln1: probably fine with what have here, but give once-over before finishing things up to determine if there needs to be any other additions (already have a few that are good things ti have here)
		DataSourcesCheckDirectionColumn <- sapply(DataSources, CheckDataSourceDirectionColumn)
		if (FALSE %in% DataSourcesCheckDirectionColumn) {
			stop(Sys.time(), " -- the following data sources have entries other than + and - in the Direction column. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckDirectionColumn], collapse=" "))
		}
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed Direction column check.", sep=""))

		DataSourcesCheckMAFIsMAF <- sapply(DataSources, CheckDataSourceMAFIsMAF)
		if (FALSE %in% DataSourcesCheckMAFIsMAF) {
			stop(Sys.time(), " -- the following data sources have variants whose MAF entry are > .5; bmass expects the MAF column to only have values <= .5. Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckMAFIsMAF], collapse=" "))
		}
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed MAF column check.", sep=""))

		DataSourcesCheckMAFFixed <- sapply(DataSources, CheckDataSourceMAFFixed) 
		if (FALSE %in% DataSourcesCheckMAFFixed) {
			stop(Sys.time(), " -- the following data sources have variants whose MAF are == 0 (or == 1); bmass expects only segregating variants (eg not fixed). Please fix and rerun bmass: ", paste(DataSources[!DataSourcesCheckMAFFixed], collapse=" "))
		} else {
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed MAF fixed check.", sep=""))
		}

		if (!is.null(ProvidedPriors)) {
			if (!is.vector(ProvidedPriors)) {
				stop(Sys.time(), " -- ProvidedPriors input is not in vector format. Please fix and rerun bmass.")
			}
			if (!is.numeric(ProvidedPriors)) {
				stop(Sys.time(), " -- ProvidedPriors input is returning false for is.numeric(). Please ensure all entries are numeric and then rerun bmass.")
			}
			if (length(ProvidedPriors) != 3^(length(DataSources)*length(SigmaAlphas))) { 
				stop(Sys.time(), " -- The number of entries in ProvidedPriors does not equal 3 ^ (the number of datasets passed to DataSources * length(SigmaAlphas)) (ie 3 ^ (", as.character(length(DataSources)), " * ", as.character(length(SigmaAlphas)), ") = ", as.character(3^(length(DataSources)*length(SigmaAlphas))), "). Please fix and rerun bmass.")
			}
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors was provided and passed checks.", sep=""))
		}
	
		#20160902 CHECK_0 -- Prob: Go through all other input variables and make sure they are the expected formats/inputs, eg numeric, character, etcetc
		if (!is.numeric(NminThreshold)) {
			stop(Sys.time(), " -- .")
		}

	} else {
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- MergedDataSources file provided, going through data checks.", sep=""))
		#20180830 NOTE -- just get rid of this little if/else section here, since not doing a 'checkmergeddatasources' thing here like I was potentially anticipating; would just do it as a separate function call, not as an if/else thing here (see commented out setup in `bmass.R` atm, and then move/alter code to potentially include something like that below/here as a placeholder for now (eg an empty function below that has that name with nothing else in it really aside from like a PH <- 1 situation or something like that))
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



