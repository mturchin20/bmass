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



#
#Annotating merged data with, if provided, GWAS SNPs
#~~~~~~

AnnotateMergedDataWithGWASsnps <- function (MergedDataSource1, GWASsnps1, BPWindow=500000) {
	GWASannot1 <- 0
	for (snpIndex in 1:nrow(GWASsnps1)) {
		if (GWASsnps1[snpIndex,]$Chr == as.numeric(as.character(MergedDataSource1["Chr"]))) {
			if (GWASsnps1[snpIndex,]$BP == as.numeric(as.character(MergedDataSource1["BP"]))) {
				GWASannot1 <- 1
			}
			else if ((GWASsnps1[snpIndex,]$BP >= as.numeric(as.character(MergedDataSource1["BP"])) - BPWindow) && (GWASsnps1[snpIndex,]$BP <= as.numeric(as.character(MergedDataSource1["BP"])) + BPWindow) && (GWASannot1 != 1)) {
				GWASannot1 <- 2	
			}
			else {
				PH <- NULL
			}
		}
			print("yaya")
			print(as.numeric(as.character(MergedDataSource1["Chr"])))
			print(GWASsnps1[snpIndex,]$Chr)
	}
	return(GWASannot1)
}
#val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE
GetZScoreAndDirection <- function(DataSources1) {
	ZScore <- qnorm(log(as.numeric(as.character(DataSources1["pValue"]))/2), lower.tail=FALSE, log.p=TRUE);
	if (DataSources1["Direction"] == "-") {
		ZScore <- ZScore * -1;
	}
	return(ZScore);
}

CheckForInfiniteZScores <- function(DataSources1_ZScores) {
	returnValue <- FALSE
	CheckValues <- is.infinite(DataSources1_ZScores)
	if (TRUE %in% CheckValues) {
		returnValue <- TRUE
	}
	return(returnValue)
}

ReplaceInfiniteZScoresWithMax <- function(DataSources1_ZScores) {
	CheckValues <- is.infinite(DataSources1_ZScores)
	if (TRUE %in% CheckValues) {
		CheckValues_positive <- is.infinite(DataSources1_ZScores) & DataSources1_ZScores > 0
		CheckValues_negative <- is.infinite(DataSources1_ZScores) & DataSources1_ZScores < 0
		maxZScore <- max(DataSources1_ZScores[!CheckValues])
		minZScore <- min(DataSources1_ZScores[!CheckValues])
		DataSources1_ZScores[CheckValues_positive] <- maxZScore
		DataSources1_ZScores[CheckValues_negative] <- minZScore
	}
	return(DataSources1_ZScores)	
}







#Conducting main bmass analyses and first-level results presentation
#~~~~~~

##collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
#20160822 CHECK_0 -- Prob: Double-check logic and go-through here
CollapseSigmaAlphasTogether <- function (inputValues1, nSigmaAlphas) { 
#	print(matrix(inputValues1, ncol=nSigmaAlphas, byrow=FALSE))
	CollapsedInputs <- apply(matrix(inputValues1, ncol=nSigmaAlphas, byrow=FALSE), 1, sum)
	return(CollapsedInputs)
}

CheckForAndReplaceOnes <- function(x) {
	returnValue1 <- x
	if (x == 1) {
		returnValue1 <- 0
	}
	return(returnValue1)
}

CheckForAndReplaceZeroes <- function(x) {
	returnValue1 <- x
	if (x == 0) {
		returnValue1 <- 1
	}
	return(returnValue1)
}

#GetMeanAcrossAlphaSigmas <- function(LogBFs1, nGammas, nSigmaAlphas) {
#        MeanAcrossAlphaSigmas <- matrix(0, ncol=ncol(LogBFs1), nrow=nGammas)
#        for (i in 1:nGammas) {
#                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
#                max <- apply(LogBFs1[SigmaAlpha_Coordinates,], 2, max)
#                LogBFs1[SigmaAlpha_Coordinates,] <- LogBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(LogBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(LogBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
#                MeanAcrossAlphaSigmas[i,] <- log10(apply(10^LogBFs1[SigmaAlpha_Coordinates,], 2, mean)) + max
#        }
#        return(MeanAcrossAlphaSigmas)
#}

####Candidate For Unit Tests####
#GetSumAcrossSigmaAlphas_withPriors(matrix(1, ncol=2, nrow=2), matrix(1, ncol=2, nrow=2), 1, 2)
#GetSumAcrossSigmaAlphas_withPriors(matrix(0, ncol=2, nrow=2), matrix(1, ncol=2, nrow=2), 1, 2)
#apply(10^matrix(0, ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)
#apply(10^matrix(c(0,1), ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)
#log10(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes))
#log10(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum))
#log10(sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes))
#matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE) - matrix(apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max), ncol=2, nrow=2, byrow=TRUE)
#log10(sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max)
#matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), ncol=2)[seq.int(1, by=2, length.out=3),]
#Test removing matrix of max values too? What about SigmaAlpa_Coordinates part?
GetSumAcrossSigmaAlphas_withPriors <- function(LogBFs1, ModelPriors, nGammas, nSigmaAlphas) {
        WeightedSumAcrossAlphaSigmas <- matrix(0, ncol=ncol(LogBFs1), nrow=nGammas)
        for (i in 1:nGammas) {
                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
                max <- apply(LogBFs1[SigmaAlpha_Coordinates,], 2, max)
                LogBFs1[SigmaAlpha_Coordinates,] <- LogBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(LogBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(LogBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
                #20160902 CHECK_0 -- Prob: Check use of max*nSigmaAlphas at end below...point is subtracting max nSigmaAlpha number of times and then summing across those rows, so shouldn't add back that max nSigmaAlpha number of times too?
		WeightedSumAcrossAlphaSigmas[i,] <- log10(sapply(apply(ModelPriors[SigmaAlpha_Coordinates,] * apply(10^LogBFs1[SigmaAlpha_Coordinates,], c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + max*nSigmaAlphas
        }
        return(WeightedSumAcrossAlphaSigmas)
}


#Convert a comma-separated string of data file locations into a vector
#20160812 20160814 CHECK_1 -- Prob: Come up with way to send off log message and then immediately leave program due to an error. Eg, here in this instance, if there are no commas the program should exit -- either the person has inputted things incorrectly or the person has supplied only a single datafile. The former should be obvious and in the latter bmass presumably cannot run on just one file/one phenotype. Soln: Found/used stop().

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
bmass <- function (DataSources, GWASsnps=NULL, SNPMarginalpValThreshold=1e-6, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), MergedDataSources=NULL, ProvidedPriors=NULL, UseFlatPriors=FALSE, PruneMarginalHits=TRUE, PruneMarginalHits_bpWindow=5e5, NminThreshold = 0, bmassSeedValue=NULL) {

	print(DataSources)

	LogFile1 <- c()
	bmassOutput <- list()
	PrevHits <- list()
	NewHits <- list()

	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning bmass.", sep=""))

	#20160823 CHECK_0: Prob -- list of Matthew functions specifically to double-check, go through, go over
	#	collapse
	#	em.priorprobs

	#Loading and checking data
	#~~~~~~

	if (is.null(MergedDataSources)) {

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning DataSources checks.", sep=""))
	
		####Candidate For Unit Tests####
		if (!is.vector(DataSources)) {
			stop(Sys.time(), " -- input variable DataSources not in vector format. bmass expects DataSources to be a vector of strings. Please fix and rerun bmass.") 
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed vector check.", sep=""))

		####Candidate For Unit Tests####
		DataSourcesCheckCharacters <- sapply(DataSources, CheckCharacterFormat)	
		if (FALSE %in% DataSourcesCheckCharacters) {
			stop(Sys.time(), " -- the following entries in DataSources were not found as characters. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckCharacters])
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed string check.", sep=""))

		####Candidate For Unit Tests####
		DataSourcesCheckExists <- sapply(DataSources, CheckVariableExists)
		if (FALSE %in% DataSourcesCheckExists) {
			stop(Sys.time(), " -- the variables associated with the following entries in DataSources were not found to exist. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckExists])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed exists check.", sep=""))
		
		#20160814 20160814 CHECK_1 -- Prob: Create way to specify which files are causing the data.frame failure here? Maybe change DataFrameCheckValues into a vector of true/false statements and then convert the trues to their text names as posible outputs? Soln: Moved to a format where returning TRUE and FALSE statements in a vector, and then pass that vector to DataSources character vector to get proper output. 
		####Candidate For Unit Tests####
		DataSourcesCheckDataFrames <- sapply(DataSources, CheckDataFrameFormat)
		if (FALSE %in% DataSourcesCheckDataFrames) {
			stop(Sys.time(), " -- the following data sources are not formatted as data.frames. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckDataFrames])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed data.frame check.", sep=""))

		####Candidate For Unit Tests####
		DataSourcesCheckHeaderNames <- CheckDataSourceHeaders(DataSources, ExpectedColumnNames)
		if (FALSE %in% DataSourcesCheckHeaderNames) {
			stop(Sys.time(), " -- the following data sources do not have all the expected column headers. The expected column headers are \"", paste(ExpectedColumnNames, collapse=" "), "\". Please fix and rerun bmass: ", DataSources[!DataSourcesCheckHeaderNames])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed column headers check.", sep=""))
		
		
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
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed Direction column check.", sep=""))

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
			LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- ProvidedPriors was provided and passed checks.", sep=""))
		}
	
		#20160902 CHECK_0 -- Prob: Go through all other input variables and make sure they are the expected formats/inputs, eg numeric, character, etcetc
		#bmass <- function (DataSources, GWASsnps=NULL, SNPMarginalpValThreshold=1e-6, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15), MergedDataSources=NULL, ProvidedPriors=NULL, UseFlatPriors=FALSE, PruneMarginalHits=TRUE, PruneMarginalHits_bpWindow=5e5, NminThreshold = 0, bmassSeedValue=NULL) {

		if (!is.numeric(NminThreshold)) {
			stop(Sys.time(), " -- .")
		}

	}
	else {
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- MergedDataSources file provided, going through data checks.", sep=""))

		#20160901 CHECK_0 -- Prob: Get this done at some point soon'ish, not soon after first rough draft completed for everything
		###### do thissssss
		#
		#Routine for checks if MergedDataSources file provided?
		#
		###### do thissssss

	}
	
	
	
	#Preparing and merging data
	#~~~~~~
	
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Beginning DataSources merging.", sep=""))

	MergedDataSources <- data.frame()
	for (CurrentDataSource in DataSources) {
	
		###### do thissssss
		#
		#Keep MAF and average across datasets? Check A1 for consistency or not expecting that since differing directions?
		#
		###### do thissssss
		
		if (nrow(MergedDataSources)==0) {
			MergedDataSources <- eval(parse(text=CurrentDataSource))	
			MergedDataSources$ZScore <- apply(MergedDataSources[,c("pValue", "Direction")], 1, GetZScoreAndDirection) 
#			MergedDataSources$pValue <- NULL
#			MergedDataSources$Direction <- NULL
			MergedDataSources_namesCurrent <- names(MergedDataSources)
			MergedDataSources_namesNew <- c()
			for (columnHeader1 in MergedDataSources_namesCurrent) {
				if (columnHeader1 %in% c("Chr", "BP", "A1", "MAF")) {
					MergedDataSources_namesNew <- c(MergedDataSources_namesNew, columnHeader1)
				}
				else {
					MergedDataSources_namesNew <- c(MergedDataSources_namesNew, paste(CurrentDataSource, "_", columnHeader1, sep=""))	
				}
			}
			names(MergedDataSources) <- MergedDataSources_namesNew
#			MergedDataSources$ChrBP <- paste(eval(parse(text=paste("MergedDataSources$", CurrentDataSource, "_Chr", sep=""))), eval(parse(text=paste("MergedDataSources$", CurrentDataSource, "_BP", sep=""))), sep="_")
			MergedDataSources$ChrBP <- paste(MergedDataSources$Chr, MergedDataSources$BP, sep="_")
		}
		else {
			CurrentDataSource_temp <- eval(parse(text=paste(CurrentDataSource, "[,c(\"Chr\", \"BP\", \"Direction\", \"pValue\", \"N\")]", sep="")))
			CurrentDataSource_temp$ZScore <- apply(CurrentDataSource_temp[,c("pValue", "Direction")], 1, GetZScoreAndDirection)
#			CurrentDataSource_temp$Direction <- NULL
#			CurrentDataSource_temp$pValue <- NULL
			CurrentDataSource_temp_namesCurrent <- names(CurrentDataSource_temp)
			CurrentDataSource_temp_namesNew <- c()
			for (columnHeader1 in CurrentDataSource_temp_namesCurrent) {
				CurrentDataSource_temp_namesNew <- c(CurrentDataSource_temp_namesNew, paste(CurrentDataSource, "_", columnHeader1, sep=""))
			}
			names(CurrentDataSource_temp) <- CurrentDataSource_temp_namesNew
			CurrentDataSource_temp$ChrBP <- paste(eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_Chr", sep=""))), eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_BP", sep=""))), sep="_")		
			eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_Chr <- NULL", sep="")))		
			eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_BP <- NULL", sep="")))		

#			print(MergedDataSources)	
#			print(CurrentDataSource_temp)
			
			###### do thissssss
			#
			# Test performances of merge function here in varying circumstances
			#
			###### do thissssss
			
			MergedDataSources <- merge(MergedDataSources, CurrentDataSource_temp, by="ChrBP")
		
			rm(CurrentDataSource_temp)
		}
	}

#	print(MergedDataSources)

	#Annotating merged data with, if provided, GWAS SNPs
	#~~~~~~

	if (is.null(GWASsnps)) {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- No GWASsnps list provided, skipping annotating MergedDataSources.", sep=""))
		MergedDataSources$GWASannot <- 0
	}
	else {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Annotating MergedDataSources with provided GWASsnps list.", sep=""))
		MergedDataSources$GWASannot <- apply(MergedDataSources, 1, AnnotateMergedDataWithGWASsnps, GWASsnps1=GWASsnps, BPWindow=500000) 
	}

	print(MergedDataSources)	

	#Calculating RSS0 and subsetting down to marginally significant SNPs
	#~~~~~~

	#Getting ZScore matrix
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Determining Z-score correlation matrix.", sep=""))
	
	ZScoresFull_CommandText <- ""
	for (DataSource in DataSources) {
		if (length(strsplit(ZScoresFull_CommandText, "")[[1]]) == 0) {
			ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, "cbind(MergedDataSources$", DataSource, "_ZScore", sep="")
		}
		else {
			ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, ",MergedDataSources$", DataSource, "_ZScore", sep="")
		}
	}
	ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, ")", sep="")
	ZScoresFull <- eval(parse(text=ZScoresFull_CommandText))

	ZScoresFullNames_CommandText <- c()
	for (DataSource in DataSources) {
		ZScoresFullNames_CommandText <- c(ZScoresFullNames_CommandText, paste(DataSource, "_ZScore", sep=""))
	}
	colnames(ZScoresFull) <- ZScoresFullNames_CommandText
	print(ZScoresFull)

	#Checking ZScore matrix for infinites and replacing with appropriate max/min values
	ZScoresFull_InfiniteCheck <- apply(ZScoresFull, 2, CheckForInfiniteZScores)
	if (TRUE %in% ZScoresFull_InfiniteCheck) {	
		warning(paste(format(Sys.time()), " -- One of your datafiles has p-values less than the threshold which R can properly convert them to log-scale, meaning their log-values return as 'Infinite'. bmass automatically replaces these 'Infinite' values with the max, non-infinite Z-scores available, but it is recommended you self-check this as well.", sep=""))
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- One of your datafiles has p-values less than the threshold which R can properly convert them to log-scale, meaning their log-values return as 'Infinite'. bmass automatically replaces these 'Infinite' values with the max, non-infinite Z-scores available, but it is recommended you self-check this as well.", sep=""))
		##UNIT TEST CANDIDATE?
		ZScoresFull <- apply(ZScoresFull, 2, ReplaceInfiniteZScoresWithMax)
		for (ZScoresFull_colName in colnames(ZScoresFull)) {
			eval(parse(text=paste("MergedDataSources$", ZScoresFull_colName, " <- ZScoresFull[,\"", ZScoresFull_colName, "\"]", sep="")))
		}
	}
#	print(ZScoresFull)

	#Getting ZScore null subset matrix, calculating correlation matrix, and associated naive multivariate statistics
	ZScoresNullSetSelection_CommandText <- ""
	for (DataSource in DataSources) {
		if (length(strsplit(ZScoresNullSetSelection_CommandText, "")[[1]]) == 0) {
			ZScoresNullSetSelection_CommandText <- paste(ZScoresNullSetSelection_CommandText, "(abs(ZScoresFull[,\"", DataSource, "_ZScore\"])<2)", sep="")
		}
		else {
			ZScoresNullSetSelection_CommandText <- paste(ZScoresNullSetSelection_CommandText, " & (abs(ZScoresFull[,\"", DataSource, "_ZScore\"])<2)", sep="") 
		}
	}
	ZScoresNullSetSelection <- eval(parse(text=ZScoresNullSetSelection_CommandText))
	ZScoresNullset <- ZScoresFull[ZScoresNullSetSelection,]
#	print(ZScoresNullset)

	ZScoresCorMatrix <- cor(ZScoresNullset)
	bmassOutput$ZScoresCorMatrix <- ZScoresCorMatrix
	
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Determining initial threshold statistics.", sep=""))
	
	ZScoresCorMatrix_Inverse <- chol2inv(chol(ZScoresCorMatrix))
	ZScoresFull_mvstat <- rowSums(ZScoresFull * (ZScoresFull %*% ZScoresCorMatrix_Inverse))
	MergedDataSources$mvstat <- ZScoresFull_mvstat
	ZScoresFull_mvstat_log10pVal <- -log10(exp(1))*pchisq(ZScoresFull_mvstat, df=ncol(ZScoresFull), log.p=TRUE, lower.tail=FALSE)
	MergedDataSources$mvstat_log10pVal <- ZScoresFull_mvstat_log10pVal

	ZScoresFull_unistat <- apply(ZScoresFull^2, 1, max)
	MergedDataSources$unistat <- ZScoresFull_unistat
	ZScoresFull_unistat_log10pVal <- -log10(exp(1))*pchisq(ZScoresFull_unistat, df=1, log.p=TRUE, lower.tail=FALSE)
	MergedDataSources$unistat_log10pVal <- ZScoresFull_unistat_log10pVal

#	print(bmassOutput)
	print(MergedDataSources)

	#Creating subset of marginally significant SNPs using SNPMarginalpValThreshold
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Subsetting down to marginally significant SNPs based on threshold: ", as.character(SNPMarginalpValThreshold) ,".", sep=""))
	MarginalHits <- MergedDataSources[MergedDataSources$mvstat_log10pVal > -log10(SNPMarginalpValThreshold) | MergedDataSources$unistat_log10pVal > -log10(SNPMarginalpValThreshold),] 
	
	print(MarginalHits)


	#Conducting main bmass analyses and first-level results presentation
	#~~~~~~

	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Conducting main bmass analysis and first-level results formatting.", sep=""))

	if (!is.null(bmassSeedValue)) {
		set.seed(bmassSeedValue)
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- setting seed with the following value: ", bmassSeedValue, ".", sep=""))	
	}
	else {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- no seed value via bmassSeedValue provided.", sep=""))	
	}
	
	ZScoresMarginal_CommandText <- ""
	for (DataSource in DataSources) {
		if (length(strsplit(ZScoresMarginal_CommandText, "")[[1]]) == 0) {
			ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, "cbind(MarginalHits$", DataSource, "_ZScore", sep="")
		}
		else {
			ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, ",MarginalHits$", DataSource, "_ZScore", sep="")
		}
	}
	ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, ")", sep="")
	ZScoresMarginal <- eval(parse(text=ZScoresMarginal_CommandText))
	
	ZScoresMarginalNames_CommandText <- c()
	for (DataSource in DataSources) {
		ZScoresMarginalNames_CommandText <- c(ZScoresMarginalNames_CommandText, paste(DataSource, "_ZScore", sep=""))
	}
	colnames(ZScoresMarginal) <- ZScoresMarginalNames_CommandText

	NsMarginal_CommandText <- ""
	for (DataSource in DataSources) {
		if (length(strsplit(NsMarginal_CommandText, "")[[1]]) == 0) {
			NsMarginal_CommandText <- paste(NsMarginal_CommandText, "cbind(MarginalHits$", DataSource, "_N", sep="")
		}
		else {
			NsMarginal_CommandText <- paste(NsMarginal_CommandText, ",MarginalHits$", DataSource, "_N", sep="")
		}
	}
	NsMarginal_CommandText <- paste(NsMarginal_CommandText, ")", sep="")
	NsMarginal <- eval(parse(text=NsMarginal_CommandText))
	NsMarginal_RowMins <- apply(NsMarginal, 1, min)
	MarginalHits$Nmin <- NsMarginal_RowMins
	
	#20160822 20160823 CHECK_1 -- Prob: Go through use of 'do.call(rbind...etc...' and double-check logic Soln: Reminder, do.call is for applying a function to a given list of arguments. Eg the contents of do.call are treated as the full set of arguments to be used, versus say an 'apply' version where the function is applied individually to each set of arguments/vectors.
	#20160822 CHECK_0 -- Prob: Change output of Matthew's code to use logBFs vs. lbf
	#20160822 CHECK_0 -- Prob: Change output of Matthew's code to match styles developed here

	print(MarginalHits)
	print(ZScoresMarginal)
	
	MarginalHits_logBFs <- compute.allBFs.fromZscores(ZScoresMarginal, ZScoresCorMatrix, MarginalHits$Nmin, MarginalHits$MAF, SigmaAlphas) 
	#MarginalHits_logBFs$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows
	MarginalHits_logBFs_Stacked <- do.call(rbind, MarginalHits_logBFs$lbf)
	
	if (!is.null(ProvidedPriors)) {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
		MarginalHits_logBFs$prior <- ProvidedPriors
	}

#	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Setting up priors.", sep=""))

#	GWASHits_EBprior <- NULL
#	FlatUnif_EBprior <- NULL
#	logBF_min	
	MarginalHits_logBFs_Stacked_AvgwPrior <- NULL
	MarginalHits_logBFs_Stacked_AvgwPrior_Min <- NULL
	Priors_Used <- NULL
	if (is.null(GWASsnps) || UseFlatPriors == TRUE) {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Setting up flat-tiered priors, GWASnps either not provided or flat prior explicitly requested.", sep=""))
	
		Prior_FlatUnif <- normalize(rep(c(0,MarginalHits_logBFs$prior[-1]),length(SigmaAlphas)))
		#nsigma=length(sigmaa)
		#origprior = rep(c(0,lbf$prior[-1]),nsigma)
		#origprior = normalize(origprior)
		
		#Prior_FlatUnif 
	
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and flat-tiered priors.", sep=""))
	
		MarginalHits_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalHits_logBFs_Stacked, Prior_FlatUnif) 
		Priors_Used <- Prior_FlatUnif	
	
		#lbf.av.origprior.glhits = lbf.av(lbf.glhits,origprior)

		#Add summary stats to marginal SNPs
		#Add SNPs x Model matrix with prior*logBFs as entries, summed (or avg'd??) across SigmaAlphas
	
	}
	else {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Setting up GWAS trained priors and analyzing GWAS hits since GWASsnps provided.", sep=""))
		GWASHits <- MarginalHits[MarginalHits$GWASannot==1,]
		GWASHits_logBFs_Stacked <- as.matrix(MarginalHits_logBFs_Stacked[,MarginalHits$GWASannot==1]) #Matrix of nSigmaAlphas x nSNPs
		
		#20160822 20160823 CHECK_1 -- Prob: Do GWAS hit analysis/work here Soln: Wrote up a few first-level GWAS hit results to get things started. Certainly will undergo further revisions down the line but fine strating point for now.
		
		Prior_GWAShitsEB <- em.priorprobs(GWASHits_logBFs_Stacked, MarginalHits_logBFs$prior, 100) #Vector with nModels*nSigmaAlphas entries
		#20160823 CHECK_0: Prob -- double check use of em.priorprobs here too with runif prior starting point
		#20160823 CHECK_0: Prob -- Do multiple runs of em.priorprobs and figure out way to compare them for consistency?
		Prior_GWAShitsEB_check2 <- em.priorprobs(GWASHits_logBFs_Stacked, MarginalHits_logBFs$prior*runif(length(MarginalHits_logBFs$prior)), 100)
	
		###### do thissssss
		#
		# Run multiple EMs to check/test for convergence? 
		#
		###### do thissssss

		Prior_GWAShitsEB_Collapsed <- CollapseSigmaAlphasTogether(Prior_GWAShitsEB, length(SigmaAlphas))
		Prior_GWAShitsEB_check2_Collapsed <- CollapseSigmaAlphasTogether(Prior_GWAShitsEB_check2, length(SigmaAlphas))

		GWASHits_PosteriorProbs <- posteriorprob(GWASHits_logBFs_Stacked, Prior_GWAShitsEB) #Matrix of nModels*nSigmaAlphas x nSNPs 
		GWASHits_PosteriorProbs_Collapsed <- apply(GWASHits_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
#		print(GWASHits_PosteriorProbs)
#		print(GWASHits_PosteriorProbs_Collapsed)
		
		PrevHits$BestModels <- apply(MarginalHits_logBFs$gamma[apply(GWASHits_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")
		GWASHits$BestModel <- apply(MarginalHits_logBFs$gamma[apply(GWASHits_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")	
#		print(PrevHits$BestModels)
		print(GWASHits)

		## this returns a list with elements  pU pD and pI
		#marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

		Prior_GWAShitsEB_ModelMatrix <- cbind(MarginalHits_logBFs$gamma, Prior_GWAShitsEB_Collapsed)[order(Prior_GWAShitsEB_Collapsed, decreasing=TRUE),]
		Prior_GWAShitsEB_ModelMatrix <- data.frame(cbind(Prior_GWAShitsEB_ModelMatrix, cumsum(Prior_GWAShitsEB_ModelMatrix[,ncol(Prior_GWAShitsEB_ModelMatrix)])))
		colnames(Prior_GWAShitsEB_ModelMatrix) <- c(DataSources, "pValue", "Cumm_pValue")
		print(Prior_GWAShitsEB_ModelMatrix)

		Prior_GWAShitsEB_ModelMatrix_AllAssoc <- apply((Prior_GWAShitsEB_ModelMatrix[,1:length(DataSources)]>0), 1, sum) == length(DataSources)
		Prior_GWAShitsEB_ModelMatrix_AllAssoc_pValSum <- sum(Prior_GWAShitsEB_ModelMatrix[Prior_GWAShitsEB_ModelMatrix_AllAssoc,ncol(Prior_GWAShitsEB_ModelMatrix)-1])
		Prior_GWAShitsEB_ModelMatrix_AllBut1Assoc <- apply((Prior_GWAShitsEB_ModelMatrix[,1:length(DataSources)]>0), 1, sum) == length(DataSources)-1
		Prior_GWAShitsEB_ModelMatrix_AllBut1Assoc_pValSum <- sum(Prior_GWAShitsEB_ModelMatrix[Prior_GWAShitsEB_ModelMatrix_AllBut1Assoc,ncol(Prior_GWAShitsEB_ModelMatrix)-1])

		Prior_GWAShitsEB_ModelMatrix_pValSupport <- c(Prior_GWAShitsEB_ModelMatrix_AllAssoc_pValSum, Prior_GWAShitsEB_ModelMatrix_AllBut1Assoc_pValSum)
		Prior_GWAShitsEB_ModelMatrix_pValSupport_Names <- c("AllAssoc", "AllBut1Assoc")
		for (DataSource in DataSources) {
			eval(parse(text=paste("Prior_GWAShitsEB_ModelMatrix_pValSupport <- c(Prior_GWAShitsEB_ModelMatrix_pValSupport, sum(Prior_GWAShitsEB_ModelMatrix[Prior_GWAShitsEB_ModelMatrix_AllBut1Assoc & Prior_GWAShitsEB_ModelMatrix[,\"", DataSource, "\"] == 0, ncol(Prior_GWAShitsEB_ModelMatrix)-1]))", sep="")))
			Prior_GWAShitsEB_ModelMatrix_pValSupport_Names <- c(Prior_GWAShitsEB_ModelMatrix_pValSupport_Names, paste("AllBut", DataSource, "Assoc", sep=""))
		}
		names(Prior_GWAShitsEB_ModelMatrix_pValSupport) <- Prior_GWAShitsEB_ModelMatrix_pValSupport_Names
	
		PrevHits$ModelCategories_CummpValues <- Prior_GWAShitsEB_ModelMatrix_pValSupport
#		print(Prior_GWAShitsEB_ModelMatrix_pValSupport)
		print(PrevHits$ModelCategories_CummpValues)
	
###?		MarginalHits_logBFs_Stacked_AvgwPrior <- ()
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and trained priors.", sep=""))
		
		MarginalHits_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalHits_logBFs_Stacked, Prior_GWAShitsEB) 
		MarginalHits_logBFs_Stacked_AvgwPrior_Min <- min(MarginalHits_logBFs_Stacked_AvgwPrior)
		Priors_Used <- Prior_GWAShitsEB
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

	NewHits <- NULL	

	print(MarginalHits_logBFs_Stacked_AvgwPrior_Min)
	
	#Determining new hits if GWASsnps were provided to determine minimum MarginalHits_logBFs_Stacked_AvgwPrior value threshold	
	#20160902 CHECK_0 -- Prob: Check that MarginalHits_logBFs_Stacked_AvgwPrior_Min is non null here?
	if (!is.null(GWASsnps)) {
		if (is.null(MarginalHits_logBFs_Stacked_AvgwPrior_Min)) {
			stop(Sys.time(), " -- MarginalHits_logBFs_Stacked_AvgwPrior_Min is NULL despite GWASsnps being provided. Unexpected error.")
		}
		NewHits <- MarginalHits[MarginalHits$GWASannot == 0 & MarginalHits$LogBFWeightedAvg >= MarginalHits_logBFs_Stacked_AvgwPrior_Min & MarginalHits$Nmin >= NminThreshold,]
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
	
	#GWASHits_PosteriorProbs <- posteriorprob(GWASHits_logBFs_Stacked, Prior_GWAShitsEB) #Matrix of nModels*nSigmaAlphas x nSNPs 
	#GWASHits_PosteriorProbs_Collapsed <- apply(GWASHits_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
	
	#MarginalHits_logBFs <- compute.allBFs.fromZscores(ZScoresMarginal, ZScoresCorMatrix, MarginalHits$Nmin, MarginalHits$MAF, SigmaAlphas) 
	##MarginalHits_logBFs$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows
	#MarginalHits_logBFs_Stacked <- do.call(rbind, MarginalHits_logBFs$lbf)

	#lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 81, 14)
	#lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
	#lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
	#lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]

#	GetSumAcrossSigmaAlphas_withPriors <- function(LogBFs1, ModelPriors, nGammas, nSigmaAlphas) {

	#Preparing final return variable bmassOutput
	bmassOutput$MarginalSNPs <- MarginalHits
	print(dim(MarginalHits_logBFs_Stacked))
	print(length(Priors_Used))
	bmassOutput$ModelPriors <- Priors_Used
	####Candidate For Unit Tests####
	print(dim(matrix(rep(Priors_Used, ncol(MarginalHits_logBFs_Stacked)), ncol=ncol(MarginalHits_logBFs_Stacked), byrow=FALSE)))
	#bmassOutput$logBFs <- GetSumAcrossSigmaAlphas_withPriors(MarginalHits_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalHits_logBFs_Stacked)), ncol=ncol(MarginalHits_logBFs_Stacked), byrow=FALSE), nrow(MarginalHits_logBFs$gamma), length(SigmaAlphas)) 
	bmassOutput$logBFs <- MarginalHits_logBFs_Stacked_SigmaAlphasSummed 
	print(dim(bmassOutput$logBFs))
	print(bmassOutput$logBFs)
	bmassOutput$PosteriorProbabilities <- MarginalHits_logBFs_Stacked_PosteriorProbabilities_Collapsed
	bmassOutput$GWASlogBFMinThreshold <- MarginalHits_logBFs_Stacked_AvgwPrior_Min
	bmassOutput$NewSNPs <- NewHits
	print(MarginalHits_logBFs_Stacked_AvgwPrior_Min)
	print(bmassOutput$NewSNPs)

#	MarginalHits_logBFs
#	MarginalHits_logBFs_Stacked

	#Returning final output variable, bmassOutput

	#return(bmassOutput)

	#lbf.av.glhits= lbf.av(lbf.glhits,ebprior.glhits)
	#nsigma=length(sigmaa)
	#origprior = rep(c(0,lbf$prior[-1]),nsigma)
	#origprior = normalize(origprior)
	#lbf.av.origprior.glhits = lbf.av(lbf.glhits,origprior)
	#lbf.uni.glhits = lbf.uni(lbf.glhits,lbf$gamma)
	#lbf.all.glhits = lbf.all(lbf.glhits,lbf$gamma)
	#
	#min(lbf.av.glhits)
	#min(lbf.av.origprior.glhits)
	##~~~
	#> min(lbf.av.glhits)
	#[1] 6.496187
	#> min(lbf.av.origprior.glhits)
	#[1] 4.723453
	#
	#lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
	#gl$lbfav = lbf.av.all
	##o = order(gl$chr, gl$pos)
	##gl = gl[o,]
	#
	#sub = gl[gl$nmin>20000,]
	#l=indephits(sub$lbfav,sub$chr,sub$pos)
	#sub=sub[l==1,]
	#newhits = sub[sub$annot==0 & sub$lbfav>5.083439 & sub$nmin>20000,c(1:4,17:28)]
	#
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
#' @param LogFile1 rbind'ing vector text output
#' @param ExpectedColumnNames Comma-separ

#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ

#if (FALSE) {
#~~~
#TimeTag 20160822
#> bmass(c("Data1", "Data2"), GWASsnps=SigSNPs)
#[1] "Data1" "Data2"
#   ChrBP Chr   BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue Data1_N
#1 1_1000   1 1000        A      0.20               -       0.0100    2500
#2 1_2000   1 2000        G      0.10               -       0.0760    2467
#3 2_3000   2 3000        T      0.06               +       0.3100    2761
#4 3_4000   3 4000        C      0.40               +       0.0056    2310
#5 4_5000   4 5000        C      0.35               -       0.7200    2632
#  Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#1   -2.5758293               -       0.0100    2500   -2.5758293         1
#2   -1.7743819               -       0.0760    2467   -1.7743819         2
#3    1.0152220               +       0.3100    2761    1.0152220         1
#4    2.7703272               +       0.0056    2310    2.7703272         0
#5   -0.3584588               -       0.7200    2632   -0.3584588         0
#     Data1_ZScore Data2_ZScore
#[1,]   -2.5758293   -2.5758293
#[2,]   -1.7743819   -1.7743819
#[3,]    1.0152220    1.0152220
#[4,]    2.7703272    2.7703272
#[5,]   -0.3584588   -0.3584588
#     Data1_ZScore Data2_ZScore
#[1,]   -1.7743819   -1.7743819
#[2,]    1.0152220    1.0152220
#[3,]   -0.3584588   -0.3584588
#             Data1_ZScore Data2_ZScore
#Data1_ZScore    1.0000000    0.5016039
#Data2_ZScore    0.5016039    1.0000000
#.
#.
#.
#  Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#1   -2.5758293               -       0.0020    2500   -3.0902323         1
#2   -1.7743819               -       0.8600    2467   -0.1763742         2
#3    1.0152220               +       0.2100    2761    1.2535654         1
#4    2.7703272               +       0.0076    2310    2.6693421         0
#5   -0.3584588               -       0.1200    2632   -1.5547736         0
#6          Inf               -       0.0000    2514         -Inf         0
#     Data1_ZScore Data2_ZScore
#[1,]   -2.5758293   -3.0902323
#[2,]   -1.7743819   -0.1763742
#[3,]    1.0152220    1.2535654
#[4,]    2.7703272    2.6693421
#[5,]   -0.3584588   -1.5547736
#[6,]          Inf         -Inf
#[1] FALSE FALSE FALSE FALSE FALSE FALSE
#[1] FALSE FALSE FALSE FALSE FALSE  TRUE
#     Data1_ZScore Data2_ZScore
#[1,]   -2.5758293   -3.0902323
#[2,]   -1.7743819   -0.1763742
#[3,]    1.0152220    1.2535654
#[4,]    2.7703272    2.6693421
#[5,]   -0.3584588   -1.5547736
#[6,]    2.7703272   -3.0902323
#             Data1_ZScore Data2_ZScore
#Data1_ZScore    1.0000000    0.5016039
#Data2_ZScore    0.5016039    1.0000000
#Warning message:
#In bmass(c("Data1", "Data2"), GWASsnps = SigSNPs) :
#  2016-08-22 13:22:29 -- One of your datafiles has p-values less than the threshold which R can properly convert them to log-scale, meaning their log-values return as 'Infinite'. bmass automatically replaces these 'Infinite' values with the max, non-infinite Z-scores available, but it is recommended you self-check this as well.
#.
#.
#.
#    ChrBP Chr    BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue Data1_N
#2  1_1000   1  1000        A      0.20               -        6e-13    2500
#4  2_3000   2  3000        T      0.06               +        3e-09    2761
#5 3_15000   3 15000        C      0.40               +        1e-07    2410
#8  4_7000   4  7000        G      0.15               +        0e+00    2514
#  Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#2    -7.200482               -        2e-13    2500    -7.348796         1
#4     5.931598               +        5e-09    2761     5.847172         1
#5     5.326724               +        4e-07    2410     5.068958         0
#8     5.931598               -        0e+00    2514    -7.348796         0
#     mvstat mvstat_log10pVal  unistat unistat_log10pVal
#2  66.29230        14.395191 54.00480         12.698970
#4  43.43998         9.432873 35.18386          8.522879
#5  33.91289         7.364091 28.37399          7.000000
#8 219.57617        47.680359 54.00480         12.698970
#     Data1_ZScore Data2_ZScore
#[1,]    -7.200482    -7.348796
#[2,]     5.931598     5.847172
#[3,]     5.326724     5.068958
#[4,]     5.931598    -7.348796
#.
#.
#.
#     Data1_ZScore Data2_ZScore
#[1,]    -7.200482    -7.348796
#[2,]     5.931598     5.847172
#[3,]     5.326724     5.068958
#[4,]     5.931598    -7.348796
#  [1]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#  [6]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [11]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [16]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [21]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [26]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [31]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [36]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [41]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [46]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [51]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [56]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [61]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [66]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [71]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [76]  0.000000e+00 2.174811e-249  0.000000e+00  0.000000e+00  0.000000e+00
# [81]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [86] 4.286876e-181  0.000000e+00  0.000000e+00 1.976263e-323  0.000000e+00
# [91]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 2.748118e-131
# [96] 1.697393e-283  0.000000e+00 7.016425e-280  0.000000e+00  0.000000e+00
#[101]  0.000000e+00  0.000000e+00  0.000000e+00  1.385787e-94 2.543654e-250
#[106]  0.000000e+00 2.648805e-245  0.000000e+00  0.000000e+00  0.000000e+00
#[111]  0.000000e+00  0.000000e+00  4.003715e-67 2.564460e-224  0.000000e+00
#[116] 8.773175e-218  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#[121]  0.000000e+00  1.000000e+00 6.643107e-155  0.000000e+00 8.553556e-143
#[126]  0.000000e+00
#      [,1] [,2]
# [1,]    0    0
# [2,]    1    0
# [3,]    2    0
# [4,]    0    1
# [5,]    1    1
# [6,]    2    1
# [7,]    0    2
# [8,]    1    2
# [9,]    2    2
#.
#.
#.
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]          [,9]         [,10]
# [1,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [2,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [3,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [4,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [5,]    0    0    0    0    0    0    0    0 4.051983e-252 3.515167e-183
# [6,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [7,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [8,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
# [9,]    0    0    0    0    0    0    0    0  0.000000e+00  0.000000e+00
#              [,11]         [,12]         [,13]         [,14]
# [1,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [2,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [3,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [4,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [5,] 7.379923e-133  9.528172e-96  5.782491e-68  1.000000e+00
# [6,] 3.501874e-286 1.149775e-252 2.169352e-226 3.099810e-156
# [7,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
# [8,] 2.021358e-282 1.714936e-247 1.085214e-219 6.197017e-144
# [9,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#
#.
#.
#.
#               [,1]          [,2]
# [1,]  0.000000e+00  0.000000e+00
# [2,]  0.000000e+00  0.000000e+00
# [3,]  0.000000e+00  0.000000e+00
# [4,]  0.000000e+00  0.000000e+00
# [5,]  1.000000e+00  1.000000e+00
# [6,] 7.404466e-157 3.099810e-156
# [7,]  0.000000e+00  0.000000e+00
# [8,] 3.359480e-145 6.197017e-144
# [9,]  0.000000e+00  0.000000e+00
#.
#.
#.
#[1] "1_1" "1_1"
#   ChrBP Chr   BP A1  MAF Data1_Direction Data1_pValue Data1_N Data1_ZScore
#2 1_1000   1 1000  A 0.20               -        6e-13    2500    -7.200482
#4 2_3000   2 3000  T 0.06               +        3e-09    2761     5.931598
#  Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot   mvstat
#2               -        2e-13    2500    -7.348796         1 66.29230
#4               +        5e-09    2761     5.847172         1 43.43998
#  mvstat_log10pVal  unistat unistat_log10pVal Nmin BestModel
#2        14.395191 54.00480         12.698970 2500       1_1
#4         9.432873 35.18386          8.522879 2761       1_1
#  Data1 Data2        pValue Cumm_pValue
#1     1     1  1.000000e+00           1
#2     1     2 8.553556e-143           1
#3     2     1 6.643107e-155           1
#4     0     0  0.000000e+00           1
#5     1     0  0.000000e+00           1
#6     2     0  0.000000e+00           1
#7     0     1  0.000000e+00           1
#8     0     2  0.000000e+00           1
#9     2     2  0.000000e+00           1
#        AllAssoc     AllBut1Assoc AllButData1Assoc AllButData2Assoc 
#               1                0                0                0 
#TimeTag 20160823 -- Program runs up to end of GWAS/prev hit analyses, moved to version 0.0.3.9000
#> GetSumAcrossSigmaAlphas_withPriors(matrix(1, ncol=2, nrow=2), matrix(1, ncol=2, nrow=2), 1, 2)
#     [,1] [,2]
#[1,]    1    1
#> GetSumAcrossSigmaAlphas_withPriors(matrix(0, ncol=2, nrow=2), matrix(1, ncol=2, nrow=2), 1, 2)
#     [,1] [,2]
#[1,]    0    0
#> apply(10^matrix(0, ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)
#     [,1] [,2]
#[1,]    0    0
#[2,]    0    0
#> 10^matrix(0, ncol=2, nrow=2)
#     [,1] [,2]
#[1,]    1    1
#[2,]    1    1
#> 10^matrix(c(0,1), ncol=2, nrow=2)
#     [,1] [,2]
#[1,]    1    1
#[2,]   10   10
#> apply(10^matrix(c(0,1), ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)
#     [,1] [,2]
#[1,]    0    0
#[2,]   10   10
#> log10(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes))
#     [,1] [,2]
#[1,] -Inf    1
#[2,] -Inf    1
#> log10(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum))
#[1]    -Inf 1.30103
#> log10(sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes))
#[1] 0.00000 1.30103
#> apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum) 
#[1]  0 20
#> sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes) 
#[1]  1 20
#> matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE) - matrix(apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max), ncol=2, nrow=2, byrow=TRUE)
#     [,1] [,2]
#[1,]    0    0
#[2,]    0    0
#> apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max)
#[1] 0 1
#> matrix(apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max), ncol=2, nrow=2, byrow=TRUE)
#     [,1] [,2]
#[1,]    0    1
#[2,]    0    1
#> log10(sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max)
#[1] 0.00000 2.30103
#> seq.int(1, by=2, length.out=3)
#[1] 1 3 5
#> matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), ncol=2)
#     [,1] [,2]
#[1,]    1    7
#[2,]    2    8
#[3,]    3    9
#[4,]    4   10
#[5,]    5   11
#[6,]    6   12
#> matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), ncol=2)[seq.int(1, by=2, length.out=3),]
#     [,1] [,2]
#[1,]    1    7
#[2,]    3    9
#[3,]    5   11
#    ChrBP Chr    BP A1  MAF Data1_Direction Data1_pValue Data1_N Data1_ZScore
#2  1_1000   1  1000  A 0.20               -        6e-13    2500    -7.200482
#4  2_3000   2  3000  T 0.06               +        3e-09    2761     5.931598
#5 3_15000   3 15000  C 0.40               +        1e-07    2410     5.326724
#8  4_7000   4  7000  G 0.15               +        0e+00    2514     5.931598
#  Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot    mvstat
#2               -        2e-13    2500    -7.348796         1  66.29230
#4               +        5e-09    2761     5.847172         1  43.43998
#5               +        4e-07    2410     5.068958         0  33.91289
#8               -        0e+00    2514    -7.348796         0 219.57617
#  mvstat_log10pVal  unistat unistat_log10pVal Nmin LogBFWeightedAvg
#2        14.395191 54.00480         12.698970 2500        12.538541
#4         9.432873 35.18386          8.522879 2761         7.411730
#5         7.364091 28.37399          7.000000 2410         5.711266
#8        47.680359 54.00480         12.698970 2514        45.346014
#       ChrBP Chr      BP A1  MAF Data1_Direction Data1_pValue Data1_N
#2     1_1000   1    1000  A 0.20               -        6e-13    2500
#4     2_3000   2    3000  T 0.06               +        3e-09    2761
#5    3_15000   3   15000  C 0.40               +        1e-07    2410
#8     4_7000   4    7000  G 0.15               +        0e+00    2514
#11 9_1038000   9 1038000  C 0.31               +        9e-23    2617
#12 9_1058000   9 1058000  C 0.29               +        5e-20    2589
#   Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#2     -7.200482               -        2e-13    2500    -7.348796         1
#4      5.931598               +        5e-09    2761     5.847172         1
#5      5.326724               +        4e-07    2410     5.068958         0
#8      9.822591               -        0e+00    2514    -7.348796         0
#11     9.822591               +        6e-23    2617     9.863372         0
#12     9.164018               +        2e-20    2589     9.262340         0
#      mvstat mvstat_log10pVal  unistat unistat_log10pVal Nmin LogBFWeightedAvg
#2   66.29230        14.395191 54.00480         12.698970 2500        12.538541
#4   43.43998         9.432873 35.18386          8.522879 2761         7.411730
#5   33.91289         7.364091 28.37399          7.000000 2410         5.711266
#8  367.95583        79.900593 96.48330         22.045757 2514        79.198156
#11 121.31709        26.343671 97.28610         22.221849 2617        24.511816
#12 106.29904        23.082543 85.79094         19.698970 2589        21.218324
#       ChrBP Chr      BP A1  MAF Data1_Direction Data1_pValue Data1_N
#2     1_1000   1    1000  A 0.20               -        6e-13    2500
#4     2_3000   2    3000  T 0.06               +        3e-09    2761
#5    3_15000   3   15000  C 0.40               +        1e-07    2410
#8     4_7000   4    7000  G 0.15               +        0e+00    2514
#11 9_1038000   9 1038000  C 0.31               +        9e-23    2617
#   Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#2     -7.200482               -        2e-13    2500    -7.348796         1
#4      5.931598               +        5e-09    2761     5.847172         1
#5      5.326724               +        4e-07    2410     5.068958         0
#8      9.822591               -        0e+00    2514    -7.348796         0
#11     9.822591               +        6e-23    2617     9.863372         0
#      mvstat mvstat_log10pVal  unistat unistat_log10pVal Nmin LogBFWeightedAvg
#2   66.29230        14.395191 54.00480         12.698970 2500        12.538541
#4   43.43998         9.432873 35.18386          8.522879 2761         7.411730
#5   33.91289         7.364091 28.37399          7.000000 2410         5.711266
#8  367.95583        79.900593 96.48330         22.045757 2514        79.198156
#11 121.31709        26.343671 97.28610         22.221849 2617        24.511816
#    ChrBP Chr    BP A1  MAF Data1_Direction Data1_pValue Data1_N Data1_ZScore
#2  1_1000   1  1000  A 0.20               -        6e-13    2500    -7.200482
#4  2_3000   2  3000  T 0.06               +        3e-09    2761     5.931598
#5 3_15000   3 15000  C 0.40               +        1e-07    2410     5.326724
#6 3_21000   3 21000  G 0.37               +        5e-08    2582     5.451310
#9  4_7000   4  7000  G 0.15               +        0e+00    2514     5.931598
#  Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot    mvstat
#2               -        2e-13    2500    -7.348796         1  66.29230
#4               +        5e-09    2761     5.847172         1  43.43998
#5               +        4e-07    2410     5.068958         0  33.91289
#6               +        9e-08    2510     5.345837         0  36.50763
#9               -        0e+00    2514    -7.348796         0 219.57617
#  mvstat_log10pVal  unistat unistat_log10pVal Nmin LogBFWeightedAvg
#2        14.395191 54.00480         12.698970 2500        12.538541
#4         9.432873 35.18386          8.522879 2761         7.411730
#5         7.364091 28.37399          7.000000 2410         5.711266
#6         7.927532 29.71679          7.301030 2510         6.257912
#9        47.680359 54.00480         12.698970 2514        45.346014
#    ChrBP Chr    BP A1  MAF Data1_Direction Data1_pValue Data1_N Data1_ZScore
#2  1_1000   1  1000  A 0.20               -        6e-13    2500    -7.200482
#4  2_3000   2  3000  T 0.06               +        3e-09    2761     5.931598
#6 3_21000   3 21000  G 0.37               +        5e-08    2582     5.451310
#9  4_7000   4  7000  G 0.15               +        0e+00    2514     5.931598
#  Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot    mvstat
#2               -        2e-13    2500    -7.348796         1  66.29230
#4               +        5e-09    2761     5.847172         1  43.43998
#6               +        9e-08    2510     5.345837         0  36.50763
#9               -        0e+00    2514    -7.348796         0 219.57617
#  mvstat_log10pVal  unistat unistat_log10pVal Nmin LogBFWeightedAvg
#2        14.395191 54.00480         12.698970 2500        12.538541
#4         9.432873 35.18386          8.522879 2761         7.411730
#6         7.927532 29.71679          7.301030 2510         6.257912
#9        47.680359 54.00480         12.698970 2514        45.346014
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


