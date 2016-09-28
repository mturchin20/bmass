context("Check DataSources")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
VectorOfTrueFalse <- c(TRUE, TRUE, TRUE, FALSE)

test_that("DataSources is a vector", {
	expect_equal(is.vector(DataSources), TRUE)
})

test_that("CheckCharacterClass checks for character class", {
	expect_equal(CheckCharacterClass("check1"), TRUE)
	expect_equal(CheckCharacterClass(12), FALSE)
	expect_equal(FALSE %in% sapply(DataSources, CheckCharacterClass), FALSE)
	expect_equal(TRUE %in% sapply(DataSources, CheckCharacterClass), TRUE)
	expect_equal(FALSE %in% VectorOfTrueFalse, TRUE)
})


#  expect_equal(str_length("a"), 1)
#  expect_equal(str_length("ab"), 2)
#  expect_equal(str_length("abc"), 3)
#
#                if (!is.vector(DataSources)) {
#                        stop(Sys.time(), " -- input variable DataSources not in vector format. bmass expects DataSources to be a vector of strings. Please fix and rerun bmass.")
#                }
#
#                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- DataSources passed vector check.", sep=""))
#
#                ####Candidate For Unit Tests####
#                DataSourcesCheckCharacters <- sapply(DataSources, CheckCharacterClass)
#                if (FALSE %in% DataSourcesCheckCharacters) {
#                        stop(Sys.time(), " -- the following entries in DataSources were not found as characters. Please fix and rerun bmass: ", DataSources[!DataSourcesCheckCharacters])
#                }
#
#CheckCharacterClass <- function (DataSource1) {
#        returnValue <- FALSE
#
#        if (is.character(DataSource1)) {
#                returnValue <- TRUE
#        }
#
#        return(returnValue)
#}





