context("Tests for Stephens2013PLoSONE.test.funcs.R")

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


rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_ZeroMatrix, bmass_TestDirection1, bmass_TestDirection2, bmass_TestDirection3, bmass_TestDirection4, envir = .GlobalEnv)

