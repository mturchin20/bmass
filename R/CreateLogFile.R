#' Checking and setting up logfile interface.
#' 
#' Description 
#' 
#' @param x Something.
#' @param y Something2.
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)

OpenLogFile <- function(FileNameBase) {
	FileNameLog <- paste(FileNameBase, ".bmass.log", sep="")
	if (file.exists(FileNameLog)) {
		file.create(FileNameLog)
		WriteLogFile(FileNameBase, paste("Warning -- file ", FileNameLog, " already existed and was overwritten.", sep=""))
	}
	else {
		file.create(FileNameLog)
	}
}

WriteLogFile <- function(FileNameBase, OutputData) {
	FileNameLog <- paste(FileNameBase, ".bmass.log", sep="")
	write(OutputData, file=FileNameLog, append=TRUE)
}

WriteTableLogFile <- function(FileNameBase, OutputData) {
	FileNameLog <- paste(FileNameBase, ".bmass.log", sep="")
	write.table(OutputData, file=FileNameLog, append=TRUE, row.names=FALSE, quote=FALSE)
}

#~~~
#file.create(..., showWarnings = TRUE)
#file.exists(...)
#file.remove(...)
#file.rename(from, to)
#file.append(file1, file2)
#file.copy(from, to, overwrite = recursive, recursive = FALSE,
#          copy.mode = TRUE, copy.date = FALSE)
#file.symlink(from, to)
#file.link(from, to)
#.
#.
#.
#file.create creates files with the given names if they do not already exist and truncates them if they do. They are created with the maximal read/write permissions allowed by the ‘umask’ setting (where relevant). By default a warning is given (with the reason) if the operation fails.
#
#file.exists returns a logical vector indicating whether the files named by its argument exist. (Here ‘exists’ is in the sense of the system's stat call: a file will be reported as existing only if you have the permissions needed by stat. Existence can also be checked by file.access, which might use different permissions and so obtain a different result. Note that the existence of a file does not imply that it is readable: for that use file.access.) What constitutes a ‘file’ is system-dependent, but should include directories. (However, directory names must not include a trailing backslash or slash on Windows.) Note that if the file is a symbolic link on a Unix-alike, the result indicates if the link points to an actual file, not just if the link exists. Lastly, note the different function exists which checks for existence of R objects.
#.
#.
#.
#~~~




