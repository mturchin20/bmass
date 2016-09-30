context("CheckDataSources")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N") 
assign("bmass_ZeroMatrix", matrix(0, nrow=2, ncol=2), envir = .GlobalEnv)
assign("bmass_TestDirection1", data.frame(Direction=c("+","-","+")), envir = .GlobalEnv)
assign("bmass_TestDirection2", data.frame(Direction=c("+",2,"+")), envir = .GlobalEnv)
assign("bmass_TestDirection3", data.frame(Direction=c("+","-","45")), envir = .GlobalEnv)
assign("bmass_TestDirection4", data.frame(Direction=c(NA,"-","+")), envir = .GlobalEnv)

test_that("CheckCharacterClass checks whether input variable is of class character", {
	expect_equal(CheckCharacterClass("check1"), TRUE)
	expect_equal(CheckCharacterClass(""), TRUE)
	expect_equal(CheckCharacterClass(12), FALSE)
	expect_equal(CheckCharacterClass(NA), FALSE)
	expect_equal(FALSE %in% sapply(DataSources, CheckCharacterClass), FALSE)
	expect_equal(TRUE %in% sapply(DataSources, CheckCharacterClass), TRUE)
})

test_that("CheckVariableExists checks whether input variable corresponds to an object that exists", {
	expect_equal(CheckVariableExists("NoDataSet"), FALSE)
	expect_equal(CheckVariableExists("bmass_ZeroMatrix"), TRUE)
	expect_equal(CheckVariableExists("bmass_TestData1"), TRUE)
	expect_equal(FALSE %in% sapply(DataSources, CheckVariableExists), FALSE)
	expect_equal(TRUE %in% sapply(DataSources, CheckVariableExists), TRUE)
})

test_that("CheckDataFrameFormat checks whether input variable corresponds to an object of datatype data.frame", {
	expect_equal(CheckDataFrameFormat("bmass_TestData1"), TRUE)
	expect_equal(CheckDataFrameFormat("bmass_ZeroMatrix"), FALSE)	
	expect_equal(FALSE %in% sapply(DataSources, CheckDataFrameFormat), FALSE)
	expect_equal(TRUE %in% sapply(DataSources, CheckDataFrameFormat), TRUE)
})

test_that("CheckDataSourceHeaders checks whether names() of object corresponding to input variables contain all expected column headers" , {
	expect_equal(ExpectedColumnNames[1] %in% names(eval(parse(text="bmass_TestData1"))), TRUE)
	expect_equal(ExpectedColumnNames[1] %in% names(eval(parse(text="bmass_ZeroMatrix"))), FALSE)
	expect_equal(TRUE %in% CheckDataSourceHeaders(DataSources, ExpectedColumnNames), TRUE)
	expect_equal(FALSE %in% CheckDataSourceHeaders(DataSources, ExpectedColumnNames), FALSE)
	expect_equal(FALSE %in% CheckDataSourceHeaders(c(DataSources, "bmass_ZeroMatrix"), ExpectedColumnNames), TRUE)
})

test_that("CheckDataSourceDirectionColumn checks whether the \'Direction\' column in the object corresponding to the input variables only contains \'+\' or \'-\'", {
	expect_equal(CheckDataSourceDirectionColumn("bmass_TestDirection1"), TRUE)
	expect_equal(CheckDataSourceDirectionColumn("bmass_TestDirection2"), FALSE)
	expect_equal(CheckDataSourceDirectionColumn("bmass_TestDirection3"), FALSE)
	expect_equal(CheckDataSourceDirectionColumn("bmass_TestDirection4"), FALSE)
	expect_equal(FALSE %in% sapply(DataSources, CheckDataSourceDirectionColumn), FALSE)
	expect_equal(FALSE %in% sapply(c(DataSources, "bmass_TestDirection1"), CheckDataSourceDirectionColumn), FALSE)
	expect_equal(TRUE %in% sapply(c(DataSources, "bmass_TestDirection2"), CheckDataSourceDirectionColumn), TRUE)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_ZeroMatrix, bmass_TestDirection1, bmass_TestDirection2, bmass_TestDirection3, bmass_TestDirection4, envir = .GlobalEnv)

