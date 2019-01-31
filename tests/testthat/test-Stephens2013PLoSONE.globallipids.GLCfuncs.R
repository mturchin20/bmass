context("Tests for Stephens2013PLoSONE.globallipids.GLCfuncs.R")

data(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs)
DataSources <- c("bmass_TestData1", "bmass_TestData2")

#test_that("indephits greedily takes top hits in a list of SNPs and removes all other SNPs within a T basepair window of them", {
#
#})

#test_that("normalize scales a vector by its sum", {
#
#})

#test_that("centered.lbf removes max from columns of log10 Bayes Factors to avoid overflow", {
#
#})

#test_that("posteriorprob calculates the posterior probability of the given models via an input of log10 Bayes Factors and the models' associated priors", {
#
#})

#test_that("em.priorprobs fits a given set of model priors by repeatedly calculating posterior probabilities over a set of SNPs and using the average posterior across SNPs as the new model priors", {
#
#})

#test_that("marginal.postprobs calculates the marginal posteriors of the multivariate categories U, D, and I given a set of posterior probabilities for all models", {
#
#})

#test_that("lbf.av calculates the average log10 Bayes Factors over a set of models weighted by those models' priors", {
#
#})

rm(bmass_TestData1, bmass_TestData2, bmass_TestSigSNPs, envir = .GlobalEnv)
