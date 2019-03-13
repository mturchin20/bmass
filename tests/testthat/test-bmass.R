context("Tests for bmass.R") 

#library("devtools"); devtools::load_all("/home/mturchin20/project/Lab_Stuff/StephensLab/bmass"); Sys.setenv(TZ='America/Chicago');
#bmass_TestData1 <- read.table("bmass_TestData1.txt", header=T)
#bmass_TestData2 <- read.table("bmass_TestData2.txt", header=T)
#bmass_TestSigSNPs <- read.table("bmass_TestSigSNPs.txt", header=T)
#load("/home/mturchin20/project/Lab_Stuff/StephensLab/bmass/data/bmass_TestData1.rda")
#load("/home/mturchin20/project/Lab_Stuff/StephensLab/bmass/data/bmass_TestData2.rda")
#load("/home/mturchin20/project/Lab_Stuff/StephensLab/bmass/data/bmass_TestSigSNPs.rda")
#DataSources <- c("bmass_TestData1", "bmass_TestData2")
#bmass_Output <- bmass(DataSources, bmass_TestSigSNPs) 
#bmass_Output2 <- bmass(DataSources, bmass_TestSigSNPs, SigmaAlphas=c(.5,.7)) 

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
bmass_TestData1 <- bmass_TestData1[bmass_TestData1$pValue > 0,]
bmass_TestData2 <- bmass_TestData2[bmass_TestData2$pValue > 0,]

assign("bmass_TestData1", bmass_TestData1, envir = .GlobalEnv)
assign("bmass_TestData2", bmass_TestData2, envir = .GlobalEnv)
assign("bmass_Output", bmass(DataSources, bmass_TestSigSNPs), envir = .GlobalEnv)
assign("bmass_Output_ZScoresCorMatrix", c(1.0000000,0.5972347,0.5972347,1.0000000), envir = .GlobalEnv)
assign("bmass_Output_GWASlogBFMinThreshold", 7.41173, envir = .GlobalEnv)
assign("bmass_Output_PreviousSNPs_SNP1", 12.53854, envir = .GlobalEnv)
assign("bmass_Output_MarginalSNPs_SNP1", 6.257912, envir = .GlobalEnv)
assign("bmass_Output2", bmass(DataSources, bmass_TestSigSNPs, SigmaAlphas=c(.5,.7)), envir = .GlobalEnv)
assign("bmass_Output2_ModelPriors_Length", 18, envir = .GlobalEnv)
assign("bmass_Output2_MarginalSNPs_logBFWAvg", c(12.218538,7.492893,5.493848), envir = .GlobalEnv)

test_that("bmass runs the main bmass function", {
	expect_equal(c(bmass_Output$ZScoresCorMatrix), bmass_Output_ZScoresCorMatrix, tolerance=1e-6)
	expect_equal(bmass_Output$GWASlogBFMinThreshold, bmass_Output_GWASlogBFMinThreshold, tolerance=1e-6)
	expect_equal(bmass_Output$PreviousSNPs$SNPs$logBFWeightedAvg[1], bmass_Output_PreviousSNPs_SNP1, tolerance=1e-6)
	expect_equal(bmass_Output$MarginalSNPs$SNPs$logBFWeightedAvg[3], bmass_Output_MarginalSNPs_SNP1, tolerance=1e-6)
	expect_equal(length(bmass_Output2$ModelPriors), bmass_Output2_ModelPriors_Length)
	expect_equal(bmass_Output2$MarginalSNPs$SNPs$logBFWeightedAvg, bmass_Output2_MarginalSNPs_logBFWAvg, tolerance=1e-6)
	expect_equal(sum(bmass_Output$MarginalSNPs$SNPs$logBFWeightedAvg) == sum(bmass_Output2$MarginalSNPs$SNPs$logBFWeightedAvg), FALSE, tolerance=1e-6)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_Output, bmass_Output_ZScoresCorMatrix, bmass_Output_GWASlogBFMinThreshold, bmass_Output_PreviousSNPs_SNP1, bmass_Output_MarginalSNPs_SNP1, bmass_Output2, bmass_Output2_ModelPriors_Length, bmass_Output2_MarginalSNPs_logBFWAvg, envir = .GlobalEnv)
