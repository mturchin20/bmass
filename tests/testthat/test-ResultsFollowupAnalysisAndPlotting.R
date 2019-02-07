context("Tests for ResultsFollowupAnalysisAndPlotting.R")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
bmass_TestData1 <- bmass_TestData1[bmass_TestData1$pValue > 0,]
bmass_TestData2 <- bmass_TestData2[bmass_TestData2$pValue > 0,]

assign("bmass_TestData1", bmass_TestData1, envir = .GlobalEnv)
assign("bmass_TestData2", bmass_TestData2, envir = .GlobalEnv)
assign("bmass_Output", bmass(DataSources, bmass_TestSigSNPs), envir = .GlobalEnv)
assign("bmass_Output_MdlPrMtrx1", c(5,8,6,1,2,3,4,7,9), envir = .GlobalEnv)
assign("bmass_Output_TopModel1", 2, envir = .GlobalEnv)
assign("bmass_Output_Marg_pI", c(7.404466e-157,3.359480e-145,3.099810e-156,6.197017e-144), envir = .GlobalEnv) 

test_that("GetModelPriorMatrix creates a matrix describing the original models and their associated trained priors", {
	expect_equal(GetModelPriorMatrix(DataSources, bmass_Output$Models, bmass_Output$ModelPriors, bmass_Output$LogFile)$ModelPriorMatrix$OrigOrder, bmass_Output_MdlPrMtrx1)	
})

test_that("GetTopModelsPerSNPViaPosteriors gets the top multivariate models per SNP based on each model's posterior probability", {
	expect_equal(as.numeric(GetTopModelsPerSNPViaPosteriors(DataSources, bmass_Output$PreviousSNPs, GetModelPriorMatrix(DataSources, bmass_Output$Models, bmass_Output$ModelPriors, bmass_Output$LogFile)$ModelPriorMatrix, bmass_Output$LogFile)$ListSNPs$TopModels[,2]), bmass_Output_TopModel1)
})

test_that("GetMarginalPosteriors gets the marginal posteriors for the U, D, and I multivariate categories for each SNP", {
	expect_equal(c(GetMarginalPosteriors(DataSources, bmass_Output$PreviousSNPs, bmass_Output$Models, bmass_Output$LogFile)$ListSNPs$Marginals$pI), bmass_Output_Marg_pI)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_Output, bmass_Output_MdlPrMtrx1, bmass_Output_TopModel1, bmass_Output_Marg_pI, envir = .GlobalEnv)
