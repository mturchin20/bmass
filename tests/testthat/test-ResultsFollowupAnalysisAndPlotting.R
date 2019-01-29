context("Tests for ResultsFollowupAnalysisAndPlotting.R")

#library("devtools"); devtools::load_all("/home/mturchin20/project/Lab_Stuff/StephensLab/bmass");
#bmass_TestData1 <- read.table("bmass_TestData1.txt", header=T)
#bmass_TestData2 <- read.table("bmass_TestData2.txt", header=T)
#bmass_TestSigSNPs <- read.table("bmass_TestSigSNPs.txt", header=T)
#DataSources <- c("bmass_TestData1", "bmass_TestData2")
#bmass_Output <- bmass(DataSources, bmass_TestSigSNPs) 

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
bmass_TestData1 <- bmass_TestData1[bmass_TestData1$pValue > 0,]
bmass_TestData2 <- bmass_TestData2[bmass_TestData2$pValue > 0,]

assign("bmass_TestData1", bmass_TestData1, envir = .GlobalEnv)
assign("bmass_TestData2", bmass_TestData2, envir = .GlobalEnv)
assign("bmass_Output", bmass(DataSources, bmass_TestSigSNPs), envir = .GlobalEnv)
#assign("bmass_Output_Marg_pI", matrix(c(7.404466e-157,3.359480e-145,3.099810e-156,6.197017e-144), ncol=2, byrow=FALSE), envir = .GlobalEnv) 
assign("bmass_Output_Marg_pI", c(7.404466e-157,3.359480e-145,3.099810e-156,6.197017e-144), envir = .GlobalEnv) 

#test_that("GetTopModelsPerSNPViaPosteriors gets the top multivariate models per SNP based on each model's posterior probability", {
#
#})

test_that("GetMarginalPosteriors gets the marginal posteriors for the U, D, and I multivariate categories for each SNP", {
	expect_equal(c(GetMarginalPosteriors(DataSources, bmass_Output$PreviousSNPs, bmass_Output$Models, bmass_Output$LogFile)$ListSNPs$Marginals$pI), bmass_Output_Marg_pI)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_Output, bmass_Output_Marg_pI, envir = .GlobalEnv)
