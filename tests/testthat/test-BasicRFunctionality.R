context("Basic R Functionality")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N") 
VectorOfTrueFalse <- c(TRUE, TRUE, FALSE)
VectorOfTrue <- c(TRUE, TRUE, TRUE)
VectorOfFalse <- c(FALSE, FALSE, FALSE)
LogTest1 <- c("Test1")
LogTest2 <- matrix(c("Test1", "Test1"), nrow=2)
assign("bmass_ZeroMatrix", matrix(0, nrow=2, ncol=2), envir = .GlobalEnv)

test_that("Check basic R functionality", {
	expect_equal(is.null(NULL), TRUE)
	expect_equal(!is.null(NULL), FALSE)
	expect_equal(is.null(NA), FALSE)
	expect_equal(is.na(NA), TRUE)
	expect_equal(is.na(435), FALSE)
	expect_equal(is.numeric(13), TRUE)
	expect_equal(is.numeric("ah"), FALSE)
	expect_equal(is.numeric(NA), FALSE)
	expect_equal(is.character("3"), TRUE)
	expect_equal(is.character(3), FALSE)
	expect_equal(is.character(NA), FALSE)
	expect_equal(is.vector(c(1,2,3)), TRUE)
	expect_equal(is.vector("2"), TRUE)
	expect_equal(is.vector(matrix(0, nrow=2, ncol=2)), FALSE)
	expect_equal(length(c(1,2,3,4,5,6,7,8,9)) != length(c(1,2)), TRUE)
	expect_equal(length(c(1,2,3,4,5,6,7,8,9)) == 3^length(c(1,2)), TRUE)
	expect_equal(length(ExpectedColumnNames) != length(DataSources), TRUE)
	expect_equal(length(ExpectedColumnNames) != 3^length(DataSources), TRUE)
	expect_equal(length(c(1,2,3,4,5,6,7,8,9)[c(1,2,3,4,5,6,7,8,9)>5]) == 4, TRUE)

	expect_error(stop("StopMessage"), "StopMessage")
	expect_error(stop(paste(LogTest1, "StopMessage", sep="")), "Test1StopMessage")
	
	expect_equal(eval(parse(text="bmass_TestData1")), bmass_TestData1)
	expect_equal(eval(parse(text=paste("bmass_TestData1", "$Direction", sep=""))), bmass_TestData1$Direction)
	expect_equal(eval(parse(text=paste("bmass_TestData1", "$Direction", sep="")))[eval(parse(text=paste("bmass_TestData1", "$Direction", sep=""))) == "+"], bmass_TestData1$Direction[bmass_TestData1$Direction == "+"])
})

test_that("DataSources is of datatype vector", {
	expect_equal(is.vector(DataSources), TRUE)
	expect_equal(!is.vector(DataSources), FALSE)
})

test_that("Check %in% functionality regarding TRUE/FALSE vectors", {
	expect_equal(TRUE %in% VectorOfTrue, TRUE)
	expect_equal(FALSE %in% VectorOfTrue, FALSE)
	expect_equal(FALSE %in% VectorOfFalse, TRUE)
	expect_equal(TRUE %in% VectorOfFalse, FALSE)
	expect_equal(TRUE %in% VectorOfTrueFalse, TRUE)
	expect_equal(FALSE %in% VectorOfTrueFalse, TRUE)
})

test_that("Check basic LogFile functionality works", {
	expect_equal(matrix(rbind(LogTest1, LogTest1)), LogTest2)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_ZeroMatrix, envir = .GlobalEnv)


