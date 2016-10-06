
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


ExploreBestModelsUsingPosteriors <- function (PreviousSNPs, MarginalSNPs, NewSNPs, Models, ModelPriors, LogFile) {

        #20160930 CHECK_0 -- Prob: Come back and go through what keeping here and what moving to later, downstream analysis section. Eg may want to keep and calculate posteriors here? AllAssoc & AllBut1Assoc may be for downstream analysis results/intrepretations?
        Prior_PreviousSNPsEB_Collapsed <- CollapseSigmaAlphasTogether(Prior_PreviousSNPsEB, length(SigmaAlphas))
        Prior_PreviousSNPsEB_check2_Collapsed <- CollapseSigmaAlphasTogether(Prior_PreviousSNPsEB_check2, length(SigmaAlphas))

        PreviousSNPs_PosteriorProbs <- posteriorprob(PreviousSNPs_logBFs_Stacked, Prior_PreviousSNPsEB) #Matrix of nModels*nSigmaAlphas x nSNPs
        PreviousSNPs_PosteriorProbs_Collapsed <- apply(PreviousSNPs_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs

        ## this returns a list with elements  pU pD and pI
        #marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

        Prior_PreviousSNPsEB_ModelMatrix <- cbind(Models, Prior_PreviousSNPsEB_Collapsed)[order(Prior_PreviousSNPsEB_Collapsed, decreasing=TRUE),]
        Prior_PreviousSNPsEB_ModelMatrix <- data.frame(cbind(Prior_PreviousSNPsEB_ModelMatrix, cumsum(Prior_PreviousSNPsEB_ModelMatrix[,ncol(Prior_PreviousSNPsEB_ModelMatrix)])))
        colnames(Prior_PreviousSNPsEB_ModelMatrix) <- c(DataSources, "pValue", "Cumm_pValue")
#               print(Prior_PreviousSNPsEB_ModelMatrix)

        Prior_PreviousSNPsEB_ModelMatrix_AllAssoc <- apply((Prior_PreviousSNPsEB_ModelMatrix[,1:length(DataSources)]>0), 1, sum) == length(DataSources)
        Prior_PreviousSNPsEB_ModelMatrix_AllAssoc_pValSum <- sum(Prior_PreviousSNPsEB_ModelMatrix[Prior_PreviousSNPsEB_ModelMatrix_AllAssoc,ncol(Prior_PreviousSNPsEB_ModelMatrix)-1])
        Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc <- apply((Prior_PreviousSNPsEB_ModelMatrix[,1:length(DataSources)]>0), 1, sum) == length(DataSources)-1
        Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc_pValSum <- sum(Prior_PreviousSNPsEB_ModelMatrix[Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc,ncol(Prior_PreviousSNPsEB_ModelMatrix)-1])

        Prior_PreviousSNPsEB_ModelMatrix_pValSupport <- c(Prior_PreviousSNPsEB_ModelMatrix_AllAssoc_pValSum, Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc_pValSum)
        Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names <- c("AllAssoc", "AllBut1Assoc")
        for (DataSource in DataSources) {
                eval(parse(text=paste("Prior_PreviousSNPsEB_ModelMatrix_pValSupport <- c(Prior_PreviousSNPsEB_ModelMatrix_pValSupport, sum(Prior_PreviousSNPsEB_ModelMatrix[Prior_PreviousSNPsEB_ModelMatrix_AllBut1Assoc & Prior_PreviousSNPsEB_ModelMatrix[,\"", DataSource, "\"] == 0, ncol(Prior_PreviousSNPsEB_ModelMatrix)-1]))", sep="")))
                Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names <- c(Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names, paste("AllBut", DataSource, "Assoc", sep=""))
        }
        names(Prior_PreviousSNPsEB_ModelMatrix_pValSupport) <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport_Names


        #20160905 NOTE -- code below eventually moved to a downstream follow-up analysis helper function. Being kept here for the moment just while focus on rewriting and cleaning up more pertinent parts of the code first
#       PrevHits <- list()
#       PrevHits$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
##      PreviousSNPs$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
##      print(PrevHits$ModelCategories_CummpValues)


        #20160905 NOTE -- Below code was not tested yet, just typing it out now to be included/tested later. This was something I wanted to start including anyways.
        #MarginalSNPs$BestModel <- apply(Models[apply(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed, 2, which.max),], 1, paste, collapse="_")
        #MarginalSNPs$BestModel_Posterior <- apply(MarginalSNPs_logBFs_Stacked_PosteriorProbabilities_Collapsed, 2, max)
        #PreviousSNPs$BestModel <- apply(Models[apply(PreviousSNPs_Posteriors_Collapsed, 2, which.max),], 1, paste, collapse="_")
        #PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_Posteriors_Collapsed, 2, max)

        ##lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 81, 14)
        ##lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
        ##lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
        ##lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]


        return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, NewSNPs=NewSNPs, LogFile=LogFile))

}




