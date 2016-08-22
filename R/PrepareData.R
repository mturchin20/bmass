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
#source("PrepareData.R")
#bmass(c("Data1", "Data2"), GWASsnps=SigSNPs)
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




#Annotating merged data with, if provided, GWAS SNPs
#~~~~~~

AnnotateMergedDataWithGWASsnps <- function (MergedDataSource1, GWASsnps1, BPWindow=500000) {
	annot1 <- 0
	for (snpIndex in 1:nrow(GWASsnps1)) {
		if (GWASsnps1[snpIndex,]$Chr == as.numeric(as.character(MergedDataSource1["Chr"]))) {
			if (GWASsnps1[snpIndex,]$BP == as.numeric(as.character(MergedDataSource1["BP"]))) {
				annot1 <- 1
			}
			else if ((GWASsnps1[snpIndex,]$BP >= as.numeric(as.character(MergedDataSource1["BP"])) - BPWindow) && (GWASsnps1[snpIndex,]$BP <= as.numeric(as.character(MergedDataSource1["BP"])) + BPWindow) && (annot1 != 1)) {
				annot1 <- 2	
			}
			else {
				PH <- NULL
			}
		}
			print("yaya")
			print(as.numeric(as.character(MergedDataSource1["Chr"])))
			print(GWASsnps1[snpIndex,]$Chr)
	}
	return(annot1)
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

#~~~
#> read.table("../data/TestData1.txt")
#   V1   V2 V3  V4        V5      V6   V7
#1 Chr   BP A1 MAF Direction pValue    N
#2   1 1000  A  .2         -     .01 2500
#3   1 2000  G  .1         -    .076 2467
#4   2 3000  T .06         +     .31 2761
#5   3 4000  C  .4         +   .0056 2310
#6   4 5000  C .35         -     .72 2632
#> Data1 <- read.table("../data/TestData1.txt")
#> CheckDataframeFormat <- function (DataFileName) {
#+         returnValue <- 0
#+
#+         if (!is.data.frame(eval(parse(text=DataFileName)))) {
#+                         returnValue <- 1
#+         }
#+
#+         return(returnValue)
#+ }
#> CheckDataframeFormat("Data1")
#[1] 0
#> Data2 <- as.matrix(Data1)
#> Data2
#     V1    V2     V3   V4    V5          V6        V7
#[1,] "Chr" "BP"   "A1" "MAF" "Direction" "pValue" "N"
#[2,] "1"   "1000" "A"  ".2"  "-"         ".01"     "2500"
#[3,] "1"   "2000" "G"  ".1"  "-"         ".076"    "2467"
#[4,] "2"   "3000" "T"  ".06" "+"         ".31"     "2761"
#[5,] "3"   "4000" "C"  ".4"  "+"         ".0056"   "2310"
#[6,] "4"   "5000" "C"  ".35" "-"         ".72"     "2632"
#> CheckDataframeFormat("Data2")
#[1] 1
#> DataList <- c("Data1", "Data2")
#> sapply(DataList, CheckDataframeFormat)
#Data1 Data2
#    0     1
#> sum(sapply(sapply(DataList, CheckDataframeFormat), as.numeric))
#[1] 1
#> sapply(DataList, CheckDataframeFormat, USE.NAMES=FALSE)
#[1] 0 1
#> sum(sapply(DataList, CheckDataframeFormat, USE.NAMES=FALSE))
#[1] 1
#
#~~~

#From http://stackoverflow.com/questions/22843775/how-to-create-periodically-send-text-to-a-log-file-while-printing-normal-outpu
#~~~
#> log_con <- file("test.log")
#> cat("write to log", file = log_con)
#> cat("write to console")
#write to console
# cat(... , file = "", sep = " ", fill = FALSE, labels = NULL,
#         append = FALSE)
#> file.create("nana2.out")
#[1] TRUE
#> FileName2 <- "nana2.out"
#> file.create(FileName2)
#[1] TRUE
#> cat("nana4\n", file = FileName2, append=TRUE)
#> cat("nana4\n", file = FileName2, append=TRUE)
#20160812 NOTE -- cat will include newlines only if I add them myself. Using write() or write.table() does automatically include newlines
#> cat("nana4\n", file = FileName3, append=TRUE)
#> write(OutputData, file=FileNameLog, append=TRUE)
#Error in cat(x, file = file, sep = c(rep.int(sep, ncolumns - 1), "\n"),  : 
#  object 'FileNameLog' not found
#> write(c(1,2,3,4), file=FileName3, append=TRUE)
#> write.table(c(1,2,3,4), file=FileName3, append=TRUE)
#Warning message:
#In write.table(c(1, 2, 3, 4), file = FileName3, append = TRUE) :
#  appending column names to file
#[  mturchin21@Michaels-MacBook-Air  ~/Documents/Work/LabMisc/StephensLab/bmass/R]$cat nana3.out 
#nana4
#nana4nana4nana4nana4
#1 2 3 4
#"x"
#"1" 1
#"2" 2
#"3" 3
#"4" 4
#~~~

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
bmass <- function (DataSources, GWASsnps=NULL, SNPMarginalpValThreshold=1e-6, ExpectedColumnNames=c("Chr", "BP", "MAF", "Direction", "pValue", "N"), SigmaAlphas = c(), MergedDataSources=NULL) {

	print(DataSources)

	LogFile1 <- c()
	bmassOutput <- list()

	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning bmass.", sep=""))

	#Loading and checking data
	#~~~~~~

	if (is.null(MergedDataSources)) {

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning DataSources checks.", sep=""))
		
		if (!is.vector(DataSources)) {
			stop(Sys.time(), " -- input variable DataSources not in vector format. bmass expects DataSources to be a vector of strings. Please fix and rerun bmass.") 
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed vector check.", sep=""))

		DataSourcesCheckCharacters <- sapply(DataSources, CheckCharacterFormat)	
		if (FALSE %in% DataSourcesCheckCharacters) {
			stop(Sys.time(), " -- the following entries in DataSources were not found as characters. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckCharacters])
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed string check.", sep=""))

		DataSourcesCheckExists <- sapply(DataSources, CheckVariableExists)
		if (FALSE %in% DataSourcesCheckExists) {
			stop(Sys.time(), " -- the variables associated with the following entries in DataSources were not found to exist. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckExists])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed exists check.", sep=""))
		
		#20160814 20160814 CHECK_1 -- Prob: Create way to specify which files are causing the data.frame failure here? Maybe change DataFrameCheckValues into a vector of true/false statements and then convert the trues to their text names as posible outputs? Soln: Moved to a format where returning TRUE and FALSE statements in a vector, and then pass that vector to DataSources character vector to get proper output. 
		DataSourcesCheckDataFrames <- sapply(DataSources, CheckDataFrameFormat)
		if (FALSE %in% DataSourcesCheckDataFrames) {
			stop(Sys.time(), " -- the following data sources are not formatted as data.frames. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckDataFrames])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed data.frame check.", sep=""))

		#Check DataSources headers for proper names

#		CheckDataSourceHeaders <- function (DataSources1, ExpectedColumnNames1) {
		DataSourcesCheckHeaderNames <- CheckDataSourceHeaders(DataSources, ExpectedColumnNames)
		if (FALSE %in% DataSourcesCheckHeaderNames) {
			stop(Sys.time(), " -- the following data sources do not have all the expected column headers. The expected column headers are \"", paste(ExpectedColumnNames, collapse=" "), "\". Please fix and rerun bmass: ", DataSources[!DataSourcesCheckHeaderNames])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed column headers check.", sep=""))
	
		DataSourcesCheckDirectionColumn <- sapply(DataSources, CheckDataSourceDirectionColumn)
		if (FALSE %in% DataSourcesCheckDirectionColumn) {
			stop(Sys.time(), " -- the following data sources have entries other than + and - in the Direction column. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckDirectionColumn])
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed Direction column check.", sep=""))

		###### do thissssss
		#
		#Check for 'X' chromosomes and convert to 23 & print warning/logfile statement
		#
		###### do thissssss

	}
	else {
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- MergedDataSources file provided, going through data checks.", sep=""))

		###### do thissssss
		#
		#Routine for checks if MergedDataSources file provided?
		#
		###### do thissssss

	}
	
	
	
	#Preparing and merging data
	#~~~~~~
	
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Beginning DataSources merging.", sep=""))

	###### do thissssss
	#
	#Change p-values based on direction of effect?
	#
	###### do thissssss
	
	#Merge different data sources into main file with all *$pValue & *$n entries for each pheno
	
	MergedDataSources <- data.frame()
	for (CurrentDataSource in DataSources) {
	
		###### do thissssss
		#
		#Keep MAF and average across datasets?
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
				#ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N")
				if (columnHeader1 %in% c("Chr", "BP")) {
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
#			CurrentDataSource_temp$ChrBP <- paste(CurrentDataSource_temp$Chr, CurrentDataSource_temp$BP, sep="_")
			CurrentDataSource_temp$ChrBP <- paste(eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_Chr", sep=""))), eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_BP", sep=""))), sep="_")		
			eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_Chr <- NULL", sep="")))		
			eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_BP <- NULL", sep="")))		

#			print(MergedDataSources)	
#			print(CurrentDataSource_temp)
			
			MergedDataSources <- merge(MergedDataSources, CurrentDataSource_temp, by="ChrBP")
		
			rm(CurrentDataSource_temp)
		}
	}

#	print(MergedDataSources)

	#Annotating merged data with, if provided, GWAS SNPs
	#~~~~~~

	if (is.null(GWASsnps)) {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Annotating MergedDataSources with provided GWASsnps list.", sep=""))
		MergedDataSources$GWASannot <- 0
	}
	else {
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- No GWASsnps list provided, skipping annotating MergedDataSources.", sep=""))
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
		ZScoresFull <- apply(ZScoresFull, 2, ReplaceInfiniteZScoresWithMax)
		for (ZScoresFull_colName in colnames(ZScoresFull)) {
			print(ZScoresFull_colName)
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
#	print(ZScoresCorMatrix)
	bmassOutput$ZScoresCorMatrix <- ZScoresCorMatrix
	
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
	MergedDataSources_MarginalHits <- MergedDataSources[MergedDataSources$mvstat_log10pVal > -log10(SNPMarginalpValThreshold) | MergedDataSources$unistat_log10pVal > -log10(SNPMarginalpValThreshold),] 
	
	print(MergedDataSources_MarginalHits)


	
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Subsetting down to marginally significant SNPs based on threshold: ", as.character(SNPMarginalpValThreshold) ,".", sep=""))


}

#Expectation
#Vector input of data.frame data sources 

#' @param DataSources input names of datasources
#' @param LogFile1 rbind'ing vector text output
#' @param ExpectedColumnNames Comma-separ

#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ

if (FALSE) {
#~~~
#> stop("nana", "nana3")
#Error: nananana3
#> stop("nana", " nana3")
#Error: nana nana3
#> val1 <- "shoop"
#> stop("nana", val1, " nana3")
#Error: nanashoop nana3
#> stop("nana ", val1, " nana3")
#Error: nana shoop nana3
#> val2 <- rbind(c("nanana"), c("nanananana"))
#> val2
#     [,1]        
#[1,] "nanana"    
#[2,] "nanananana"
#> stop("nana ", val1, " nana3 ", val2)
#Error: nana shoop nana3 nananananananana
#> print(val2)
#     [,1]        
#[1,] "nanana"    
#[2,] "nanananana"
#> time()
#Error in time.default() : argument "x" is missing, with no default
#> sys.time()
#Error: could not find function "sys.time"
#> Sys.time()
#[1] "2016-08-14 01:11:17 CDT"
#> Sys.time()
#[1] "2016-08-14 01:11:18 CDT"
#> paste(Sys.time(), " nanana"))
#Error: unexpected ')' in "paste(Sys.time(), " nanana"))"
#> paste(Sys.time(), " nanana")
#[1] "2016-08-14 01:11:34  nanana"
#> paste(Sys.time(), "nanana")
#[1] "2016-08-14 01:11:39 nanana"
#> paste(Sys.time(), " -- nanana")
#[1] "2016-08-14 01:11:44  -- nanana"
#> rbind(paste(Sys.time(), " -- nanana"), paste(Sys.time(), " -- nanana2"))
#     [,1]                             
#[1,] "2016-08-14 01:11:59  -- nanana" 
#[2,] "2016-08-14 01:11:59  -- nanana2"
#> rbind(paste(Sys.time(), " -- nanana"), paste(Sys.time(), " -- nanana2"))
#     [,1]                             
#[1,] "2016-08-14 01:12:04  -- nanana" 
#[2,] "2016-08-14 01:12:04  -- nanana2"
#> stop(Sys.time())
#Error: 2016-08-14 01:18:16
#> stop(Sys.time(), " -- blahlbah")
#Error: 2016-08-14 01:18:21 -- blahlbah
#> myfunc <- function(v1) {
#+   deparse(substitute(v1))
#+ }
#> deparse(substitute(val1))
#[1] "val1"
#> substitute(val1) 
#val1
#> val1
#[1] "shoop"
#> val1 <- c(1,2,3,4,5)
#> val1[c(TRUE,FALSE,FALSE,TRUE,FALSE)]
#[1] 1 4
#> TRUE %in% c(FALSE, FALSE, FALSE)
#[1] FALSE
#> TRUE %in% c(FALSE, FALSE, TRUE)
#[1] TRUE
#> FALSE %in% c(TRUE, TRUE, TRUE)
#[1] FALSE
#> FALSE %in% c(TRUE, TRUE, FALSE)
#[1] TRUE
#> blah <- data.frame()
#> blah
#data frame with 0 columns and 0 rows
#> nrow(blah)
#[1] 0
#> Data1
#  Chr   BP A1  MAF Direction pValue    N
#1   1 1000  A 0.20         -  0.0100 2500
#2   1 2000  G 0.10         -  0.0760 2467
#3   2 3000  T 0.06         +  0.3100 2761
#4   3 4000  C 0.40         +  0.0056 2310
#5   4 5000  C 0.35         -  0.7200 2632
#> blah1 <- data.frame()
#> blah1 <- Data1
#> blah1
#  Chr   BP A1  MAF Direction pValue    N
#1   1 1000  A 0.20         -  0.0100 2500
#2   1 2000  G 0.10         -  0.0760 2467
#3   2 3000  T 0.06         +  0.3100 2761
#4   3 4000  C 0.40         +  0.0056 2310
#5   4 5000  C 0.35         -  0.7200 2632
#> eval(parse(text="blah1$nana1 <- c(1,2,3,4,5)"))
#> blah1
#  Chr   BP A1  MAF Direction pValue    N nana1
#1   1 1000  A 0.20         -  0.0100 2500     1
#2   1 2000  G 0.10         -  0.0760 2467     2
#3   2 3000  T 0.06         +  0.3100 2761     3
#4   3 4000  C 0.40         +  0.0056 2310     4
#5   4 5000  C 0.35         -  0.7200 2632     5
#> eval(parse(text=paste("blah1", "$nana2 <- c(1,2,3,4,5)", sep="")))
#> blah1
#  Chr   BP A1  MAF Direction pValue    N nana1 nana2
#1   1 1000  A 0.20         -  0.0100 2500     1     1
#2   1 2000  G 0.10         -  0.0760 2467     2     2
#3   2 3000  T 0.06         +  0.3100 2761     3     3
#4   3 4000  C 0.40         +  0.0056 2310     4     4
#5   4 5000  C 0.35         -  0.7200 2632     5     5
#> val2 <- "shana"
#> paste(val2, " haha ", paste("blahblah", val2), sep="")
#[1] "shana haha blahblah shana"
#> paste(val2, " haha ", paste("blahblah", val2, sep="_"), sep="")
#[1] "shana haha blahblah_shana"
#> blah1$ChrBP <- paste(blah1$Chr, blah1$BP, sep="_")
#> blah1
#  Chr   BP A1  MAF Direction pValue    N nana1 nana2  ChrBP
#1   1 1000  A 0.20         -  0.0100 2500     1     1 1_1000
#2   1 2000  G 0.10         -  0.0760 2467     2     2 1_2000
#3   2 3000  T 0.06         +  0.3100 2761     3     3 2_3000
#4   3 4000  C 0.40         +  0.0056 2310     4     4 3_4000
#5   4 5000  C 0.35         -  0.7200 2632     5     5 4_5000
#~~~
#20160815 -- program runs up to merge
#> source("PrepareData.R")
#> bmass(c("Data1", "Data2"))
#[1] "Data1" "Data2"
#   ChrBP Data1_Chr Data1_BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue
#1 1_1000         1     1000        A      0.20               -       0.0100
#2 1_2000         1     2000        G      0.10               -       0.0760
#3 2_3000         2     3000        T      0.06               +       0.3100
#4 3_4000         3     4000        C      0.40               +       0.0056
#5 4_5000         4     5000        C      0.35               -       0.7200
#  Data1_N Data2_pValue Data2_N
#1    2500       0.0100    2500
#2    2467       0.0760    2467
#3    2761       0.3100    2761
#4    2310       0.0056    2310
#5    2632       0.7200    2632
#> bmass(c("Data1", "Data2", "Data3"))
#[1] "Data1" "Data2" "Data3"
#   ChrBP Data1_Chr Data1_BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue
#1 1_1000         1     1000        A      0.20               -       0.0100
#2 1_2000         1     2000        G      0.10               -       0.0760
#3 2_3000         2     3000        T      0.06               +       0.3100
#4 3_4000         3     4000        C      0.40               +       0.0056
#5 4_5000         4     5000        C      0.35               -       0.7200
#  Data1_N Data2_pValue Data2_N Data3_pValue Data3_N
#1    2500       0.0100    2500       0.0100    2500
#2    2467       0.0760    2467       0.0760    2467
#3    2761       0.3100    2761       0.3100    2761
#4    2310       0.0056    2310       0.0056    2310
#5    2632       0.7200    2632       0.7200    2632
#
#> bmass(c("Data1", "Data2"))
#[1] "Data1" "Data2"
#   ChrBP Chr   BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue Data1_N
#1 1_1000   1 1000        A      0.20               -       0.0100    2500
#2 1_2000   1 2000        G      0.10               -       0.0760    2467
#3 2_3000   2 3000        T      0.06               +       0.3100    2761
#4 3_4000   3 4000        C      0.40               +       0.0056    2310
#5 4_5000   4 5000        C      0.35               -       0.7200    2632
#  Data2_pValue Data2_N
#1       0.0100    2500
#2       0.0760    2467
#3       0.3100    2761
#4       0.0056    2310
#5       0.7200    2632
#> SigSNPs
#  Chr   BP
#1   1 1000
#2   2 3000
#> bmass(c("Data1", "Data2"), GWASsnps=SigSNPs)
#[1] "Data1" "Data2"
#   ChrBP Chr   BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue Data1_N
#1 1_1000   1 1000        A      0.20               -       0.0100    2500
#2 1_2000   1 2000        G      0.10               -       0.0760    2467
#3 2_3000   2 3000        T      0.06               +       0.3100    2761
#4 3_4000   3 4000        C      0.40               +       0.0056    2310
#5 4_5000   4 5000        C      0.35               -       0.7200    2632
#  Data2_pValue Data2_N GWASannot
#1       0.0100    2500         1
#2       0.0760    2467         2
#3       0.3100    2761         1
#4       0.0056    2310         0
#5       0.7200    2632         0
#> val6
#list()
#> val6$bps <- c(1,2,3,4,5)
#> val6
#$bps
#[1] 1 2 3 4 5
#
#> val6$yal1 <- matrix(c(1,2,3,4), nrow=2) 
#> val6
#$bps
#[1] 1 2 3 4 5
#
#$yal1
#     [,1] [,2]
#[1,]    1    3
#[2,]    2    4
#
#> bmass(c("Data1", "Data2"), GWASsnps=SigSNPs)
#[1] "Data1" "Data2"
#   ChrBP Chr   BP Data1_A1 Data1_MAF Data1_N Data1_ZScore Data2_N Data2_ZScore
#1 1_1000   1 1000        A      0.20    2500    2.3263479    2500    2.3263479
#2 1_2000   1 2000        G      0.10    2467    1.4325027    2467    1.4325027
#3 2_3000   2 3000        T      0.06    2761   -0.4958503    2761   -0.4958503
#4 3_4000   3 4000        C      0.40    2310   -2.5363960    2310   -2.5363960
#5 4_5000   4 5000        C      0.35    2632   -0.5828415    2632   -0.5828415
#  GWASannot
#1         1
#2         2
#3         1
#4         0
#5         0
#> bmass(c("Data1", "Data2"), GWASsnps=SigSNPs)
#[1] "Data1" "Data2"
#   ChrBP Chr   BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue Data1_N
#1 1_1000   1 1000        A      0.20               -       0.0100    2500
#2 1_2000   1 2000        G      0.10               -       0.0760    2467
#3 2_3000   2 3000        T      0.06               +       0.3100    2761
#4 3_4000   3 4000        C      0.40               +       0.0056    2310
#5 4_5000   4 5000        C      0.35               -       0.7200    2632
#  Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#1    2.3263479               -       0.0100    2500    2.3263479         1
#2    1.4325027               -       0.0760    2467    1.4325027         2
#3   -0.4958503               +       0.3100    2761   -0.4958503         1
#4   -2.5363960               +       0.0056    2310   -2.5363960         0
#5   -0.5828415               -       0.7200    2632   -0.5828415         0
#> bmass(c("Data1", "Data2"), GWASsnps=SigSNPs)
#[1] "Data1" "Data2"
#   ChrBP Chr   BP Data1_A1 Data1_MAF Data1_Direction Data1_pValue Data1_N
#1 1_1000   1 1000        A      0.20               -       0.0100    2500
#2 1_2000   1 2000        G      0.10               -       0.0760    2467
#3 2_3000   2 3000        T      0.06               +       0.3100    2761
#4 3_4000   3 4000        C      0.40               +       0.0056    2310
#5 4_5000   4 5000        C      0.35               -       0.7200    2632
#  Data1_ZScore Data2_Direction Data2_pValue Data2_N Data2_ZScore GWASannot
#1   -2.3263479               -       0.0100    2500   -2.3263479         1
#2   -1.4325027               -       0.0760    2467   -1.4325027         2
#3    0.4958503               +       0.3100    2761    0.4958503         1
#4    2.5363960               +       0.0056    2310    2.5363960         0
#5    0.5828415               -       0.7200    2632    0.5828415         0
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
}


