#' Checking and preparing input datafiles.
#' 
#' Description 
#' 
#' @param x Something.
#' @param y Something2.
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)

#Data1 <- read.table("../data/TestData1.txt")

ExpectedColumnNames1 <- c("Chr", "BP", "A1", "MAF", "Direction", "p_Value", "N")

CheckDataframeFormat <- function (DataFileName, logFile) {
	returnValue <- 0
	
	if (!is.data.frame(eval(parse(text=DataFileName)))) {
			returnValue <- 1
	}

	return(returnValue)
}

~~~
> read.table("../data/TestData1.txt")
   V1   V2 V3  V4        V5      V6   V7
1 Chr   BP A1 MAF Direction p_Value    N
2   1 1000  A  .2         -     .01 2500
3   1 2000  G  .1         -    .076 2467
4   2 3000  T .06         +     .31 2761
5   3 4000  C  .4         +   .0056 2310
6   4 5000  C .35         -     .72 2632
> Data1 <- read.table("../data/TestData1.txt")
> CheckDataframeFormat <- function (DataFileName) {
+         returnValue <- 0
+
+         if (!is.data.frame(eval(parse(text=DataFileName)))) {
+                         returnValue <- 1
+         }
+
+         return(returnValue)
+ }
> CheckDataframeFormat("Data1")
[1] 0
> Data2 <- as.matrix(Data1)
> Data2
     V1    V2     V3   V4    V5          V6        V7
[1,] "Chr" "BP"   "A1" "MAF" "Direction" "p_Value" "N"
[2,] "1"   "1000" "A"  ".2"  "-"         ".01"     "2500"
[3,] "1"   "2000" "G"  ".1"  "-"         ".076"    "2467"
[4,] "2"   "3000" "T"  ".06" "+"         ".31"     "2761"
[5,] "3"   "4000" "C"  ".4"  "+"         ".0056"   "2310"
[6,] "4"   "5000" "C"  ".35" "-"         ".72"     "2632"
> CheckDataframeFormat("Data2")
[1] 1
> DataList <- c("Data1", "Data2")
> sapply(DataList, CheckDataframeFormat)
Data1 Data2
    0     1
> sum(sapply(sapply(DataList, CheckDataframeFormat), as.numeric))
[1] 1
> sapply(DataList, CheckDataframeFormat, USE.NAMES=FALSE)
[1] 0 1
> sum(sapply(DataList, CheckDataframeFormat, USE.NAMES=FALSE))
[1] 1

~~~

#From http://stackoverflow.com/questions/22843775/how-to-create-periodically-send-text-to-a-log-file-while-printing-normal-outpu
~~~
> log_con <- file("test.log")
> cat("write to log", file = log_con)
> cat("write to console")
write to console
 cat(... , file = "", sep = " ", fill = FALSE, labels = NULL,
         append = FALSE)
~~~

CheckColumns <- function (ExpectedColumnNames, DataFile, logFile) {
	DataColumnNames <- names(DataFile)
	for ColumnName in ExpectedColumnNames {
		if (! ColumnName %in% DataColumnNames) {

		}
	}
}

PrepareData <- function (ColumnNames, DataFileNames, logFile) {

	DataFrameCheckValues <- apply(DataFileNames, CheckDataframeFormate)
	if (sum(apply(DataFrameCheckValues, as.numeric)) > 0) {
		cat("
	}

}


