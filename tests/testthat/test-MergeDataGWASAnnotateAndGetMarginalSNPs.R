context("Tests for MergeDataGWASAnnotateAndGetMarginalSNPs.R")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
ExpectedColumnNames <- c("Chr", "BP", "A1", "MAF", "Direction", "pValue", "N") 

assign("AnnotateTest1", c(1,2,1,0,0,0,0,0,0,0,0), envir = .GlobalEnv)
assign("AnnotateTest2", c(1,0,1,0,0,0,0,0,0,0,0), envir = .GlobalEnv)
assign("ZScoreOutput1", c(-7.2004822,-1.7743819,5.9315982,2.7703272,5.3267239,5.4513104,-0.3584588,Inf,0.5100735,-2.1570727,1.6448536), envir = .GlobalEnv)
assign("ZScoreOutput2", c(-7.3487961,-0.1763742,5.8471721,2.6693421,5.0689577,5.3458374,-1.5547736,-Inf,0.4399132,-1.9514798,1.1263911), envir = .GlobalEnv)
assign("ZScoreOutput1Replace", c(-7.2004822,-1.7743819,5.9315982,2.7703272,5.3267239,5.4513104,-0.3584588,5.9315982,0.5100735,-2.1570727,1.6448536), envir = .GlobalEnv)
assign("ZScoreOutput2Replace", c(-7.3487961,-0.1763742,5.8471721,2.6693421,5.0689577,5.3458374,-1.5547736,-7.3487961,0.4399132,-1.9514798,1.1263911), envir = .GlobalEnv)

test_that("AnnotateDataWithGWASSNPs annotates input list of SNPs with a second list of GWAS SNPs", {
	expect_equal(AnnotateDataWithGWASSNPs(bmass_TestData1, bmass_TestSigSNPs), AnnotateTest1)
	expect_equal(AnnotateDataWithGWASSNPs(bmass_TestData1, bmass_TestSigSNPs, BPWindow=1), AnnotateTest2)
})

test_that("GetZScoreAndDirection extracts p-value and direction information and converts to ZScore per SNP", {
	expect_equal(apply(bmass_TestData1, 1, GetZScoreAndDirection), ZScoreOutput1)
	expect_equal(apply(bmass_TestData2, 1, GetZScoreAndDirection), ZScoreOutput2)
})

test_that("CheckForInfiniteZScores checks whether there is an infinite value among the input", {
	expect_equal(TRUE %in% sapply(ZScoreOutput1, CheckForInfiniteZScores), TRUE)
	expect_equal(TRUE %in% sapply(ZScoreOutput2[1:5], CheckForInfiniteZScores), FALSE)
})

test_that("ReplaceInfiniteZScoresWithMax replaces infinite ZScores with the maximum or minimum (as appropriate re: direction of effect) present among other ZScores", {
	expect_equal(ReplaceInfiniteZScoresWithMax(ZScoreOutput1), ZScoreOutput1Replace)
	expect_equal(ReplaceInfiniteZScoresWithMax(ZScoreOutput2), ZScoreOutput2Replace)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, AnnotateTest1, AnnotateTest2, ZScoreOutput1, ZScoreOutput2, ZScoreOutput1Replace, ZScoreOutput2Replace, envir = .GlobalEnv)
