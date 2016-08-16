#' Checking and preparing input datafiles.
#' 
#' Description 
#' 
#' @param x Something.
#' @param y Something2.
####' @param DataFileList Comma-separated list of all the datafile names being analyzed (these are assumed to be the names of the phenotypes being analyzed). This should be in the same order as DataFileLocations. The default value for this paramter is NULL. If NULL, the list of datafile names are derived from the list of datafile locations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList. The default value for this parameter is NULL. Users must supply a list of datafile locations to run bmass.
#' @param ExpectedColumnNames Comma-separated list of the expected column names to be found at the head of each datafile. The default value for this parameter is Chr,BP,A1,MAF,Direction,p_Value,N. Users should not supply or alter this parameter.
#' @param OutputFileBase
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)

Data1 <- read.table("../data/TestData1.txt")
Data1
ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "p_Value", "N")
Data2 <- as.matrix(Data1)
Data2
DataList <- c("Data1", "Data2")
#	if (!is.data.frame(eval(parse(text=DataFileName)))) {

CheckCharacterFormat <- function (DataSource1) {
	returnValue <- FALSE
	
	if (!is.character(DataSource1)) {
		returnValue <- TRUE
	}

	return(returnValue)
}

CheckVariableExists <- function (DataSource1) {
	returnValue <- FALSE
	
	if (!exists(DataSource1)) {
		returnValue <- TRUE
	}

	return(returnValue)
}

CheckDataFrameFormat <- function (DataSource1) {
	returnValue <- FALSE
	
	if (!is.data.frame(DataSource1)) {
			returnValue <- TRUE
	}

	return(returnValue)
}

CheckDataSourceHeaders <- function (DataSources1, ExpectedColumnNames1) {
	returnValueVector <- c()
	for DataSource1 in DataSources1 {
		returnValue <- TRUE
		DataSource1ColumnNames <- names(eval(parse(text=DataSource1)))
		for ColumnName in ExpectedColumnNames {
			if (! ColumnName %in% DataColumnNames) {
				returnValue <- FALSE
			}
		}
		returnValueVector <- c(returnValueVector, returnValue)
	}

	return(returnValueVector)
}


#~~~
#> read.table("../data/TestData1.txt")
#   V1   V2 V3  V4        V5      V6   V7
#1 Chr   BP A1 MAF Direction p_Value    N
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
#[1,] "Chr" "BP"   "A1" "MAF" "Direction" "p_Value" "N"
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
#20160812 CHECK_0 -- Prob: Come up with way to send off log message and then immediately leave program due to an error. Eg, here in this instance, if there are no commas the program should exit -- either the person has inputted things incorrectly or the person has supplied only a single datafile. The former should be obvious and in the latter bmass presumably cannot run on just one file/one phenotype.

~~~
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
~~~

~~~
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
~~~



#This is going to be the main function that goes through each of the steps from beginning to end. Hypothetically, all the other functions presented here should be used through the PrepareData process (or as a subfunction of one of the functions being used in PrepareData)
#PrepareData <- function (ColumnNames, DataFileNames, OutputFileBase) {
#PrepareData <- function (ExpectedColumnNames, DataFileNames, DataFileLocations, OutputFileBase) { #20160812 NOTE -- deciding to remove 'DataFileNames' as a requested/given input. Just going to assume that we are given a list of file locations and I will parse that list for the 'names' of the datasources 
#PrepareData <- function (ExpectedColumnNames, DataFileLocations, OutputFileBase) { #20160814 NOTE -- Changing direction and just assuming input is a single vector that contains all the proper data.frame datasources and working from there. Final output will be a list that has all the output. 'Logfile' will just be a variable included in final list output, developed by continula 'rbind' calls with text output additions. Also deciding to move 'PrepareData' to just a 'MainWorkFlow' or 'Main' that I'll dev in each sub R package and then eventually move to a main source.  
bmass <- function (DataSources, ExpectedColumnNames=c(c("Chr", "BP", "MAF", "Direction", "p_Value", "N"), MergedDataSources=NULL) {

	LogFile1 <- c()

	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning bmass."))

	#Loading and checking data
	#~~~~~~

	if (is.null(MergedDataSources)) {

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- beginning DataSources checks."))
		
		if (!is.vector(DataSources)) {
			stop(Sys.time(), " -- input variable DataSources not in vector format. bmass expects DataSources to be a vector of strings. Please fix and rerun bmass.") 
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed vector check."))

		DataSourcesCheckCharacters <- sapply(DataSources, CheckCharacterFormat)	
		if (FALSE %in% DataSourcesCheckCharacters) {
			stop(Sys.time(), " -- the following entries in DataSources were not found as characters. Please fix and rerun bmass: ", DataSources[DataSourcesCheckCharacters])
		}
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed string check."))

		DataSourcesCheckExists <- sapply(DataSources, CheckVariableExists)
		if (FALSE %in% DataSourcesCheckExists) {
			stop(Sys.time(), " -- the variables associated with the following entries in DataSources were not found to exist. Please fix and rerun bmass: ", DataSources[DataSourcesCheckExists])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed exists check."))
		
		eval(parse(text=DataFileName))

		#20160814 20160814 CHECK_1 -- Prob: Create way to specify which files are causing the data.frame failure here? Maybe change DataFrameCheckValues into a vector of true/false statements and then convert the trues to their text names as posible outputs? Soln: Moved to a format where returning TRUE and FALSE statements in a vector, and then pass that vector to DataSources character vector to get proper output. 
		DataSourcesCheckDataFrames <- sapply(DataSources, CheckDataFrameFormat)
		if (FALSE %in% DataSourcesCheckDataFrames) {
			stop(Sys.time(), " -- the following data sources are not formatted as data.frames. Please fix and rerun bmass: ", DataSources[DataSourcesCheckDataFrames])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed data.frame check."))

		#Check DataSources headers for proper names

		CheckDataSourceHeaders <- function (DataSources1, ExpectedColumnNames1) {
		DataSourcesCheckHeaderNames <- CheckDataSourceHeaders(DataSources, ExpectedColumnNames)
		if (FALSE %in% DataSourcesCheckHeaderNames) {
			stop(Sys.time(), " -- the following data sources do not have all the expected column headers. The expected column headers are ", paste(ExpectedColumnNames, sep=" "), ". Please fix and rerun bmass: ", DataSources[DataSourcesCheckHeaderNames])
		}

		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- DataSources passed column headers check."))

	}
	else {
		
		LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- MergedDataSources file provided, going through data checks."))

		###### do thissssss
		#
		#Routine for checks if MergedDataSources file provided?
		#
		###### do thissssss

	}

	#Preparing and merging data
	~~~~~~
	
	LogFile1 <- rbind(LogFile1, paste(format(Sys.time()), " -- Beginning DataSources merging."))

	###### do thissssss
	#
	#Change p-values based on direction of effect?
	#
	###### do thissssss

	#Merge different data sources into main file with all *$pValue & *$n entries for each pheno
	
	MergedDataSources <- data.frame()
	for CurrentDataSource in DataSources {
		
		eval(parse(text=DataSource1))	
		ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "p_Value", "N")
	
		if (nrow(MergedDataSources)==0) {
			MergedDataSources <- eval(parse(text=CurrentDataSource))	
		}
		else {
		
			CurrentDataSource_temp <- eval(parse(text=paste(CurrentDataSource, "[,c(\"))
			CurrentDataSource_temp$ChrBP <- paste(CurrentDataSource_temp$Ch	

			eval(parse(text=paste(CurrentDataSource, "$ChrBP <- ", paste(paste(CurrentDataSource, "$Chr", sep=""), paste(CurrentDataSource, "$BP", sep=""), sep="_"), sep="")))
			> blah1$ChrBP <- paste(blah1$Chr, blah1$BP, sep="_")
			
			merge(MergedDataSources, eval(parse(text=CurrentDataSource))[,c("",)], by="rsID")
		}
	
		rm(CurrentDataSource_temp)
	
	}





}

#Expectation
#Vector input of data.frame data sources 

#' @param DataSources input names of datasources
#' @param LogFile1 rbind'ing vector text output
#' @param ExpectedColumnNames Comma-separ

#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ


~~~
> stop("nana", "nana3")
Error: nananana3
> stop("nana", " nana3")
Error: nana nana3
> val1 <- "shoop"
> stop("nana", val1, " nana3")
Error: nanashoop nana3
> stop("nana ", val1, " nana3")
Error: nana shoop nana3
> val2 <- rbind(c("nanana"), c("nanananana"))
> val2
     [,1]        
[1,] "nanana"    
[2,] "nanananana"
> stop("nana ", val1, " nana3 ", val2)
Error: nana shoop nana3 nananananananana
> print(val2)
     [,1]        
[1,] "nanana"    
[2,] "nanananana"
> time()
Error in time.default() : argument "x" is missing, with no default
> sys.time()
Error: could not find function "sys.time"
> Sys.time()
[1] "2016-08-14 01:11:17 CDT"
> Sys.time()
[1] "2016-08-14 01:11:18 CDT"
> paste(Sys.time(), " nanana"))
Error: unexpected ')' in "paste(Sys.time(), " nanana"))"
> paste(Sys.time(), " nanana")
[1] "2016-08-14 01:11:34  nanana"
> paste(Sys.time(), "nanana")
[1] "2016-08-14 01:11:39 nanana"
> paste(Sys.time(), " -- nanana")
[1] "2016-08-14 01:11:44  -- nanana"
> rbind(paste(Sys.time(), " -- nanana"), paste(Sys.time(), " -- nanana2"))
     [,1]                             
[1,] "2016-08-14 01:11:59  -- nanana" 
[2,] "2016-08-14 01:11:59  -- nanana2"
> rbind(paste(Sys.time(), " -- nanana"), paste(Sys.time(), " -- nanana2"))
     [,1]                             
[1,] "2016-08-14 01:12:04  -- nanana" 
[2,] "2016-08-14 01:12:04  -- nanana2"
> stop(Sys.time())
Error: 2016-08-14 01:18:16
> stop(Sys.time(), " -- blahlbah")
Error: 2016-08-14 01:18:21 -- blahlbah
> myfunc <- function(v1) {
+   deparse(substitute(v1))
+ }
> deparse(substitute(val1))
[1] "val1"
> substitute(val1) 
val1
> val1
[1] "shoop"
> val1 <- c(1,2,3,4,5)
> val1[c(TRUE,FALSE,FALSE,TRUE,FALSE)]
[1] 1 4
> TRUE %in% c(FALSE, FALSE, FALSE)
[1] FALSE
> TRUE %in% c(FALSE, FALSE, TRUE)
[1] TRUE
> FALSE %in% c(TRUE, TRUE, TRUE)
[1] FALSE
> FALSE %in% c(TRUE, TRUE, FALSE)
[1] TRUE
> blah <- data.frame()
> blah
data frame with 0 columns and 0 rows
> nrow(blah)
[1] 0
> Data1
  Chr   BP A1  MAF Direction p_Value    N
1   1 1000  A 0.20         -  0.0100 2500
2   1 2000  G 0.10         -  0.0760 2467
3   2 3000  T 0.06         +  0.3100 2761
4   3 4000  C 0.40         +  0.0056 2310
5   4 5000  C 0.35         -  0.7200 2632
> blah1 <- data.frame()
> blah1 <- Data1
> blah1
  Chr   BP A1  MAF Direction p_Value    N
1   1 1000  A 0.20         -  0.0100 2500
2   1 2000  G 0.10         -  0.0760 2467
3   2 3000  T 0.06         +  0.3100 2761
4   3 4000  C 0.40         +  0.0056 2310
5   4 5000  C 0.35         -  0.7200 2632
> eval(parse(text="blah1$nana1 <- c(1,2,3,4,5)"))
> blah1
  Chr   BP A1  MAF Direction p_Value    N nana1
1   1 1000  A 0.20         -  0.0100 2500     1
2   1 2000  G 0.10         -  0.0760 2467     2
3   2 3000  T 0.06         +  0.3100 2761     3
4   3 4000  C 0.40         +  0.0056 2310     4
5   4 5000  C 0.35         -  0.7200 2632     5
> eval(parse(text=paste("blah1", "$nana2 <- c(1,2,3,4,5)", sep="")))
> blah1
  Chr   BP A1  MAF Direction p_Value    N nana1 nana2
1   1 1000  A 0.20         -  0.0100 2500     1     1
2   1 2000  G 0.10         -  0.0760 2467     2     2
3   2 3000  T 0.06         +  0.3100 2761     3     3
4   3 4000  C 0.40         +  0.0056 2310     4     4
5   4 5000  C 0.35         -  0.7200 2632     5     5
> val2 <- "shana"
> paste(val2, " haha ", paste("blahblah", val2), sep="")
[1] "shana haha blahblah shana"
> paste(val2, " haha ", paste("blahblah", val2, sep="_"), sep="")
[1] "shana haha blahblah_shana"
> blah1$ChrBP <- paste(blah1$Chr, blah1$BP, sep="_")
> blah1
  Chr   BP A1  MAF Direction p_Value    N nana1 nana2  ChrBP
1   1 1000  A 0.20         -  0.0100 2500     1     1 1_1000
2   1 2000  G 0.10         -  0.0760 2467     2     2 1_2000
3   2 3000  T 0.06         +  0.3100 2761     3     3 2_3000
4   3 4000  C 0.40         +  0.0056 2310     4     4 3_4000
5   4 5000  C 0.35         -  0.7200 2632     5     5 4_5000
~~~





