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
	
		if (!is.numeric(NminThreshold)) {
			stop(Sys.time(), " -- input variable NminThreshold is not numeric as needed. Please fix and rerub bmass.")
		}

	} else {
		
		LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- MergedDataSources file provided, going through data checks.", sep=""))
	
	}

	return(LogFile)
}
