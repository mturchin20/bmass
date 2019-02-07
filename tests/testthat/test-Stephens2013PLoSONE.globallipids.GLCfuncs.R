context("Tests for Stephens2013PLoSONE.globallipids.GLCfuncs.R")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")
bmass_TestData1 <- bmass_TestData1[bmass_TestData1$pValue > 0,]
bmass_TestData2 <- bmass_TestData2[bmass_TestData2$pValue > 0,]

assign("bmass_TestData1", bmass_TestData1, envir = .GlobalEnv)
assign("bmass_TestData2", bmass_TestData2, envir = .GlobalEnv)
assign("bmass_Output", bmass(DataSources, bmass_TestSigSNPs), envir = .GlobalEnv)
bmass_Output_Posteriors <- bmass_Output$PreviousSNPs$Posteriors[,3:4]; bmass_Output_Posteriors[5,1] <- .45; bmass_Output_Posteriors[6,1] <- .2; bmass_Output_Posteriors[8,1] <- .35; bmass_Output_Posteriors[5,2] <- .45; bmass_Output_Posteriors[6,2] <- .35; bmass_Output_Posteriors[8,2] <- .2;
assign("bmass_Output_Posteriors", bmass_Output_Posteriors, envir = .GlobalEnv)
assign("IndepTest1",  matrix(c(1e-10,1e-7,1e-5,1e-3,.9,1,1,1,1,1,250000,200000,1,650000,25100), ncol=3, byrow=FALSE), envir = .GlobalEnv)
assign("IndepTest1_Output1", c(1,0,0,0,0), envir = .GlobalEnv)
assign("IndepTest1_Output2", c(1,0,1,1,0), envir = .GlobalEnv)
assign("IndepTest1_Normalize1", c(2.222023e-01,1.777618e-01,8.888091e-07,5.777259e-01,2.230911e-02), envir = .GlobalEnv)
assign("IndepTest1_CenteredLBF", c(0,-3,-5,-7,-9.954243), envir = .GlobalEnv)
assign("IndepTest1_PosteriorProb", c(9.979737e-01,1.995947e-03,2.993921e-05,3.991895e-07,5.544298e-10), envir = .GlobalEnv)
assign("IndepTest1_EM.PriorProbs", c(1,2e-15,3e-25,4e-35,8.467544e-50), envir = .GlobalEnv)
assign("bmass_Output_Posteriors_pIMarg", c(0.20,0.35), envir = .GlobalEnv)
assign("IndepTest1_lbfAvg", 10.30191, envir = .GlobalEnv)

test_that("indephits greedily takes top hits in a list of SNPs and removes all other SNPs within a T basepair window of them", {
	expect_equal(indephits(-log10(IndepTest1[,1]), IndepTest1[,2], IndepTest1[,3]), IndepTest1_Output1)
	expect_equal(indephits(-log10(IndepTest1[,1]), IndepTest1[,2], IndepTest1[,3], T=100000), IndepTest1_Output2)
})

test_that("normalize scales a vector by its sum", {
	expect_equal(normalize(IndepTest1[,3]), IndepTest1_Normalize1, tolerance = 1e-6)
})

test_that("centered.lbf removes max from columns of log10 Bayes Factors to avoid overflow", {
	expect_equal(centered.lbf(matrix(rep(-log10(IndepTest1[,1]), 3), ncol=3, byrow=FALSE))[,1], IndepTest1_CenteredLBF, tolerance = 1e-6)
})

test_that("posteriorprob calculates the posterior probability of the given models via an input of log10 Bayes Factors and the models' associated priors", {
	expect_equal(posteriorprob(matrix(rep(-log10(IndepTest1[,1]), 3), ncol=3, byrow=FALSE), c(2,4,6,8,10))[,1], IndepTest1_PosteriorProb, tolerance = 1e-6)
})

test_that("em.priorprobs fits a given set of model priors by repeatedly calculating posterior probabilities over a set of SNPs and using the average posterior across SNPs as the new model priors", {
	expect_equal(em.priorprobs(matrix(rep(-log10(IndepTest1[,1]), 3), ncol=3, byrow=FALSE), c(2,4,6,8,10), n=5), IndepTest1_EM.PriorProbs, tolerance = 1e-6)
})

test_that("marginal.postprobs calculates the marginal posteriors of the multivariate categories U, D, and I given a set of posterior probabilities for all models", {
	expect_equal(marginal.postprobs(bmass_Output_Posteriors, bmass_Output$Models, 1)$pI[,1], bmass_Output_Posteriors_pIMarg, tolerance = 1e-6)
})

test_that("lbf.av calculates the average log10 Bayes Factors over a set of models weighted by those models' priors", {
	expect_equal(lbf.av(matrix(rep(-log10(IndepTest1[,1]), 3), ncol=3, byrow=FALSE), c(2,4,6,8,10))[1], IndepTest1_lbfAvg, tolerance = 1e-6)
})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, bmass_Output_Posteriors, IndepTest1, IndepTest1_Output1, IndepTest1_Output2, IndepTest1_Normalize1, IndepTest1_CenteredLBF, IndepTest1_PosteriorProb, IndepTest1_EM.PriorProbs, bmass_Output_Posteriors_pIMarg, IndepTest1_lbfAvg, envir = .GlobalEnv)
