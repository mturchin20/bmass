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
ExpectedColumnNames1 <- c("Chr", "BP", "A1", "MAF", "Direction", "p_Value", "N")
Data2 <- as.matrix(Data1)
Data2
DataList <- c("Data1", "Data2")


CheckDataframeFormat <- function (DataFileName, OutputFileBase) {
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
VectorizeDataLocations <- function(DataFileLocations) {

	VectorizedDataFileLocations <- c()

	if (grep(",", DataFileLocations)) {

	}
	else {
		log file ( explain error then end with 'exiting...' )
		exit program? or supply a flag?
		problem -- at this stage we are trying to create the name for the log file so have not set it up yet
		use warning() and then exit program? say something like 'cannot initialize function/program'
		create log file or just create a variable that contains log message output and supply it as one of the output variables via list()? eg output$log
		can just keep doing rbind() rbind()....start each log message also with some sort of R time call?
	}

	return(VectorizedDataFileLocations)

}

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

#20160812 NOTE -- Check/force 'character' setting of input here?
GetListOfFileNames <- function(DataFileLocations) {

	

}

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

#I don't recall if I have a greater function for 'OpenData' other than
OpenData <- function (DataFileNames, DataFileLocations, OutputFileBase) {
	
	#Check whether OutputFileBase is NULL and if so, create it from 
	OutputFileBase <- 
		
	#Get list of 	

	return(OutputFileBase)

}

CheckColumns <- function (ExpectedColumnNames, DataFile, OutputFileBase) {
	DataColumnNames <- names(DataFile)
	for ColumnName in ExpectedColumnNames {
		if (! ColumnName %in% DataColumnNames) {

		}
	}
}

#This is going to be the main function that goes through each of the steps from beginning to end. Hypothetically, all the other functions presented here should be used through the PrepareData process (or as a subfunction of one of the functions being used in PrepareData)
#PrepareData <- function (ColumnNames, DataFileNames, OutputFileBase) {
#PrepareData <- function (ExpectedColumnNames, DataFileNames, DataFileLocations, OutputFileBase) { #20160812 NOTE -- deciding to remove 'DataFileNames' as a requested/given input. Just going to assume that we are given a list of file locations and I will parse that list for the 'names' of the datasources 
PrepareData <- function (ExpectedColumnNames, DataFileLocations, OutputFileBase) {

	#First thing should be setup log file
	SetUpLogFile <-

	DataFileLocations <- VectorizeDataLocations(DataFileLocations)

	DataFileNames <- GetListOfFilesNames(DataFileLocations)

	OutputFileBase <- OpenData(DataFileNames, DataFileLocations, OutputFileBase) 
	#If no outputfilebase is provided figure it out from the first input file given?

	DataFrameCheckValues <- apply(DataFileNames, CheckDataframeFormat)
	if (sum(apply(DataFrameCheckValues, as.numeric)) > 0) {
		WriteLogFile(OutputFileBase, paste("	
		cat("
	
	}

}





#' @param DataFileNames Comma-separated list of all the phenotypes (ie names associated with each datafile) being analyzed. This should be in the same order as DataFileLocations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList.
#' @param ExpectedColumnNames Comma-separ

