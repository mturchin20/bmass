#' Checking and preparing input datafiles.
#' 
#' Description 
#' 
#' @param x Something.
#' @param y Something2.
####' @param DataFileList Comma-separated list of all the datafile names being analyzed (these are assumed to be the names of the phenotypes being analyzed). This should be in the same order as DataFileLocations. The default value for this paramter is NULL. If NULL, the list of datafile names are derived from the list of datafile locations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList. The default value for this parameter is NULL. Users must supply a list of datafile locations to run bmass.
#' @param ExpectedColumnNames Comma-separated list of the expected column names to be found at the head of each datafile. The default value for this parameter is Chr,BP,A1,MAF,Direction,p_Value,N. Users should not supply or alter this parameter.
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)

Data1 <- read.table("../data/TestData1.txt")
Data1
ExpectedColumnNames1 <- c("Chr", "BP", "A1", "MAF", "Direction", "p_Value", "N")
Data2 <- as.matrix(Data1)
Data2
DataList <- c("Data1", "Data2")


CheckDataframeFormat <- function (DataFileName, FileNameBase) {
	returnValue <- 0
	
	if (!is.data.frame(eval(parse(text=DataFileName)))) {
			returnValue <- 1
	}

	return(returnValue)
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
#~~~

OpenData <- function (DataFileNames, DataFileLocations) {
	FileNameBase
	#Get list of 	

}

CheckColumns <- function (ExpectedColumnNames, DataFile, FileNameBase) {
	DataColumnNames <- names(DataFile)
	for ColumnName in ExpectedColumnNames {
		if (! ColumnName %in% DataColumnNames) {

		}
	}
}

#This is going to be the main function that goes through each of the steps from beginning to end. Hypothetically, all the other functions presented here should be used through the PrepareData process (or as a subfunction of one of the functions being used in PrepareData)
#PrepareData <- function (ColumnNames, DataFileNames, FileNameBase) {
PrepareData <- function (ExpectedColumnNames, DataFileNames, DataFileLocations) {

	FileNameBase <- OpenData(DataFileNames, DataFileLocations) 

	DataFrameCheckValues <- apply(DataFileNames, CheckDataframeFormat)
	if (sum(apply(DataFrameCheckValues, as.numeric)) > 0) {
		WriteLogFile(FileNameBase, paste("	
		cat("
	
	}

}





#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ

