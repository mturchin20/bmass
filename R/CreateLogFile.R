OpenLogFile <- function(OutputFileBase) {
	FileNameLog <- paste(OutputFileBase, ".bmass.log", sep="")
	if (file.exists(FileNameLog)) {
		file.create(FileNameLog)
		WriteLogFile(OutputFileBase, paste("Warning -- file ", FileNameLog, " already existed and was overwritten.", sep=""))
	}
	else {
		file.create(FileNameLog)
	}
}

WriteLogFile <- function(OutputFileBase, OutputData) {
	FileNameLog <- paste(OutputFileBase, ".bmass.log", sep="")
	write(OutputData, file=FileNameLog, append=TRUE)
}

#' WriteTableLogFile
#'
#' WriteTableLogFile
#' 
#' @keywords internal
#'
#' @importFrom utils write.table
WriteTableLogFile <- function(OutputFileBase, OutputData) {
	FileNameLog <- paste(OutputFileBase, ".bmass.log", sep="")
	write.table(OutputData, file=FileNameLog, append=TRUE, row.names=FALSE, quote=FALSE)
}
