context("Tests for GetResultsFromMarginalSNPsAndFormat.R")

set.seed(8345)
assign("CollapseTest1", c(12,15,18), envir = .GlobalEnv)
assign("CheckOnesTest1", c(0,2,3,4,5,6,7,8,9), envir = .GlobalEnv)
assign("CheckOnesTest2", c(0,0,0,0), envir = .GlobalEnv)
assign("CheckZeroesTest1", c(1,2,3,4,5,6,7,8,9), envir = .GlobalEnv)
assign("CheckZeroesTest2", c(0,10,0,10), envir = .GlobalEnv)
assign("ModelPriors", rep(c(.5,.35,.15), 3), envir = .GlobalEnv)
assign("MarginalSNPs_logBFs_Stacked", matrix(replicate(3, sample(c(1,2,3,4,5,6,7,8,9))), nrow=9, byrow=TRUE), envir = .GlobalEnv)
assign("CollapseTest_wPriors", c(6.000000,6.544502,2.176091,2.698970,0.544068,4.217484,1.000000,1.845098,2.217484), envir = .GlobalEnv)
assign("CollapseTest_wPriors_Basic1", c(1,3,5,7,9,11), envir = .GlobalEnv)

test_that("CollapseSigmaAlphasTogether sums multiple entries of the same 'model' over different sigma_alpha hyperparameter values", {
	expect_equal(CollapseSigmaAlphasTogether(matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow=FALSE), 3), CollapseTest1)
})

test_that("CheckForAndReplaceOnes checks for 1 entries and replaces them with 0s", {
	expect_equal(c(apply(matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow=FALSE), c(1,2), CheckForAndReplaceOnes)), CheckOnesTest1)
	expect_equal(c(apply(10^matrix(0, ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)), CheckOnesTest2) 
})

test_that("CheckForAndReplaceZeroes checks for 0 entries and replaces them with 1s", {
	expect_equal(c(apply(apply(matrix(c(1,2,3,4,5,6,7,8,9), ncol=3, byrow=FALSE), c(1,2), CheckForAndReplaceOnes), c(1,2), CheckForAndReplaceZeroes)), CheckZeroesTest1)
	expect_equal(c(apply(10^matrix(c(0,1), ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)), CheckZeroesTest2) 
})

test_that("GetSumAcrossSigmaAlphas_withPriors sums multiple entries (in log10) of the same 'model' over different sigma_alpha hyperparameter values while multiplying each entry by that model's prior", {
	expect_equal(c(GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(ModelPriors, ncol(MarginalSNPs_logBFs_Stacked)), nrow=length(ModelPriors), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), 3, 3)), CollapseTest_wPriors, tolerance=1e-6)
	expect_equal(c(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), ncol=2)[seq.int(1, by=2, length.out=3),]), CollapseTest_wPriors_Basic1) 
})

rm(CollapseTest1, CheckOnesTest1, CheckOnesTest2, CheckZeroesTest1, CheckZeroesTest2, ModelPriors, MarginalSNPs_logBFs_Stacked, CollapseTest_wPriors, CollapseTest_wPriors_Basic1, envir = .GlobalEnv)
