#' Checking and preparing input datafiles.
#' 
#' Description 
#' 
#' @param x Something.
#' @param y Something2.
####' @param DataFileList Comma-separated list of all the datafile names being analyzed (these are assumed to be the names of the phenotypes being analyzed). This should be in the same order as DataFileLocations. The default value for this paramter is NULL. If NULL, the list of datafile names are derived from the list of datafile locations.
#' @param DataFileLocations Comma-separated list of all the datafile locations being analyzed. This should be in the same order as PhenotypeList. The default value for this parameter is NULL. Users must supply a list of datafile locations to run bmass.
#' @param ExpectedColumnNames Comma-separated list of the expected column names to be found at the head of each datafile. The default value for this parameter is Chr,BP,A1,MAF,Direction,p_Value,N. Users should not supply or alter this parameter.
#' @param OutputFileBase
#' @return A merged and combined dataset The sum of \code{x} and \code{y}.
#' @examples
#' func(1, 1)
#' func(10, 1)

FinalizeAndFormatResults <- function(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile) {

	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and trained priors.", sep=""))

        #Pruning marginal hits by LogBFWeightedAvg if requested
        if (PruneMarginalSNPs == TRUE) {
                #20160901 CHECK_0 -- Prob: Go over indepthits function, rewrite, or just lightly edit? redo names, double-check functionality? def get some unit testing in there
                #MarginalSNPs <- MarginalSNPs[indephits(MarginalSNPs$LogBFWeightedAvg, MarginalSNPs$Chr, MarginalSNPs$BP, T=PruneMarginalSNPs_bpWindow)==1,]
                MarginalSNPs_PrunedList <- indephits(MarginalSNPs$LogBFWeightedAvg, MarginalSNPs$Chr, MarginalSNPs$BP, T=PruneMarginalSNPs_bpWindow)
                MarginalSNPs <- MarginalSNPs[MarginalSNPs_PrunedList==1,]
                MarginalSNPs_logBFs$lbf <- lapply(MarginalSNPs_logBFs$lbf, function(x) { return(x[,MarginalSNPs_PrunedList==1]) })
                MarginalSNPs_logBFs_Stacked <- MarginalSNPs_logBFs_Stacked[,MarginalSNPs_PrunedList==1]
        }

        print(MarginalSNPs)
#       print(summary(MarginalSNPs_logBFs$lbf))

        #Preparing log Bayes factors matrix, Gammas x SNPs
        #bmassOutput$logBFs <- GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalSNPs_logBFs_Stacked)), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), nrow(MarginalSNPs_logBFs$gamma), length(SigmaAlphas))
        MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed <- GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalSNPs_logBFs_Stacked)), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), nrow(MarginalSNPs_logBFs$gamma), length(SigmaAlphas))
        MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed <- cbind(MarginalSNPs_logBFs$gamma, MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed)
        colnames(MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed) <- c(DataSources, MarginalSNPs$ChrBP)

        print(MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed)
        print(Priors_Used)
        print(MarginalSNPs_logBFs_Stacked)

        #Preparing posterior probabilities, Gammas x SNPs
        MarginalSNPs_logBFs_Stacked_PosteriorProbabilities <- posteriorprob(MarginalSNPs_logBFs_Stacked, Priors_Used)
        MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed <- apply(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas))
        MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed <- cbind(MarginalSNPs_logBFs$gamma, MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed)
        colnames(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed) <- c(DataSources, MarginalSNPs$ChrBP)
        #MarginalSNPs_PosteriorProbabilities <- apply(posteriorprob(MarginalSNPs_logBFs_Stacked, Priors_Used), 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas))

        print(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed)

        #20160905 CHECK_0 -- Prob: Convert either 'MarginalSNPs_logBFs_Stacked_PosteriorProbabilities' to 'PosteriorProbs' or 'PreviousSNPs_PosteriorProbs' to 'PreviousSNPs_PosteriorProbabilities'
        #20160905 CHECK_0 -- Prob: Change 'MarginalSNPs' to 'MarginalSNPs' probably?
        #20160905 NOTE -- Below code was not tested yet, just typing it out now to be included/tested later. This was something I wanted to start including anyways.
        #MarginalSNPs$BestModel <- apply(MarginalSNPs_logBFs$gamma[apply(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed, 2, which.max),], 1, paste, collapse="_")
        #MarginalSNPs$BestModel_Posterior <- apply(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed, 2, max)

        #PreviousSNPs_PosteriorProbs <- posteriorprob(PreviousSNPs_logBFs_Stacked, Prior_PreviousSNPsEB) #Matrix of nModels*nSigmaAlphas x nSNPs
        #PreviousSNPs_PosteriorProbs_Collapsed <- apply(PreviousSNPs_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
        #PreviousSNPs$BestModel <- apply(MarginalSNPs_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")
        #PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, max)

        #lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 81, 14)
        #lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
        #lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
        #lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]

        #20160905 CHECK_0 -- Prob: Move this below intialization section to proper beginning of .R file code as necessary once change/reorganization occurs
        NewSNPs <- NULL

        print(MarginalSNPs_logBFs_Stacked_AvgwPrior_Min)

        #Determining new hits if GWASsnps were provided to determine minimum MarginalSNPs_logBFs_Stacked_AvgwPrior value threshold
        #20160902 CHECK_0 -- Prob: Check that MarginalSNPs_logBFs_Stacked_AvgwPrior_Min is non null here?
        if (!is.null(GWASsnps)) {
                if (is.null(MarginalSNPs_logBFs_Stacked_AvgwPrior_Min)) {
                        stop(Sys.time(), " -- MarginalSNPs_logBFs_Stacked_AvgwPrior_Min is NULL despite GWASsnps being provided. Unexpected error.")
                }
                NewSNPs <- MarginalSNPs[MarginalSNPs$GWASannot == 0 & MarginalSNPs$LogBFWeightedAvg >= MarginalSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$Nmin >= NminThreshold,]
        }



        #Preparing final return variable bmassOutput
        bmassOutput$MarginalSNPs$SNPs <- MarginalSNPs
        print(dim(MarginalSNPs_logBFs_Stacked))
        print(length(Priors_Used))
        bmassOutput$ModelPriors <- Priors_Used
        ####Candidate For Unit Tests####
        print(dim(matrix(rep(Priors_Used, ncol(MarginalSNPs_logBFs_Stacked)), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE)))
        #bmassOutput$logBFs <- GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(Priors_Used, ncol(MarginalSNPs_logBFs_Stacked)), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), nrow(MarginalSNPs_logBFs$gamma), length(SigmaAlphas))
        bmassOutput$MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed
        print(dim(bmassOutput$logBFs))
        print(bmassOutput$logBFs)
        bmassOutput$MarginalSNPs$Posteriors <- MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed
        bmassOutput$GWASlogBFMinThreshold <- MarginalSNPs_logBFs_Stacked_AvgwPrior_Min
        bmassOutput$NewSNPs$SNPs <- NewSNPs
        print(MarginalSNPs_logBFs_Stacked_AvgwPrior_Min)
        print(bmassOutput$NewSNPs)

        #PreviousSNPs$BestModel <- apply(MarginalSNPs_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")
        #PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, max)

        #bmassOutput <- list()
        #bmassOutput$ModelPriors <- NULL
        #bmassOutput$MarginalSNPs <- list()
        #bmassOutput$MarginalSNPs$logBFs <- NULL
        #bmassOutput$MarginalSNPs$Posteriors <- NULL
        #bmassOutput$NewSNPs <- list()
        #bmassOutput$NewSNPs$SNPs <- NULL
        #bmassOutput$NewSNPs$logBFs <- NULL
        #bmassOutput$NewSNPs$Posteriors <- NULL
        #bmassOutput$GWASlogBFMinThreshold <- NULL
        #bmassOutput$PreviousSNPs <- list()
        #bmassOutput$PreviousSNPs$SNPs <- NULL
        #bmassOutput$PreviousSNPs$logBFs <- NULL
        #bmassOutput$PreviousSNPs$Posteriors <- NULL


        #Returning final output variable, bmassOutput

        #return(bmassOutput)

        ##extract lbfs for all the new hits
        #lbf.newhits= lbf.bigmat[,gl$nmin>20000]
        #lbf.newhits= lbf.newhits[,l==1]
        #lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>5.083439 & sub$nmin>20000]
        #pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
        #pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))
        #
        #ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
        #pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
        #bestclass= apply(pp.newhits.classmatrix, 2,which.max)
        #cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))







}















