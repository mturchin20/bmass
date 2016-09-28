
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

#20160812 NOTE -- main file for collecting functions that operate on organized, finished bmass output

#Data1 <- read.table("../data/TestData1.txt")
#Data1
#ExpectedColumnNames1 <- c("Chr", "BP", "A1", "MAF", "Direction", "p_Value", "N")
#Data2 <- as.matrix(Data1)
#Data2

#Conducting main bmass analyses and first-level results presentation
#~~~~~~

##collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
#20160822 CHECK_0 -- Prob: Double-check logic and go-through here
CollapseSigmaAlphasTogether <- function (inputValues1, nSigmaAlphas) {
#       print(matrix(inputValues1, ncol=nSigmaAlphas, byrow=FALSE))
        CollapsedInputs <- apply(matrix(inputValues1, ncol=nSigmaAlphas, byrow=FALSE), 1, sum)
        return(CollapsedInputs)
}

CheckForAndReplaceOnes <- function(x) {
        returnValue1 <- x
        if (x == 1) {
                returnValue1 <- 0
        }
        return(returnValue1)
}

CheckForAndReplaceZeroes <- function(x) {
        returnValue1 <- x
        if (x == 0) {
                returnValue1 <- 1
        }
        return(returnValue1)
}

#GetMeanAcrossAlphaSigmas <- function(LogBFs1, nGammas, nSigmaAlphas) {
#        MeanAcrossAlphaSigmas <- matrix(0, ncol=ncol(LogBFs1), nrow=nGammas)
#        for (i in 1:nGammas) {
#                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
#                max <- apply(LogBFs1[SigmaAlpha_Coordinates,], 2, max)
#                LogBFs1[SigmaAlpha_Coordinates,] <- LogBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(LogBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(LogBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
#                MeanAcrossAlphaSigmas[i,] <- log10(apply(10^LogBFs1[SigmaAlpha_Coordinates,], 2, mean)) + max
#        }
#        return(MeanAcrossAlphaSigmas)
#}

####Candidate For Unit Tests####
#GetSumAcrossSigmaAlphas_withPriors(matrix(1, ncol=2, nrow=2), matrix(1, ncol=2, nrow=2), 1, 2)
#GetSumAcrossSigmaAlphas_withPriors(matrix(0, ncol=2, nrow=2), matrix(1, ncol=2, nrow=2), 1, 2)
#apply(10^matrix(0, ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)
#apply(10^matrix(c(0,1), ncol=2, nrow=2), c(1,2), CheckForAndReplaceOnes)
#log10(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes))
#log10(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum))
#log10(sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes))
#matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE) - matrix(apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max), ncol=2, nrow=2, byrow=TRUE)
#log10(sapply(apply(apply(10^matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + apply(matrix(c(0,1), ncol=2, nrow=2, byrow=TRUE), 2, max)
#matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), ncol=2)[seq.int(1, by=2, length.out=3),]
#Test removing matrix of max values too? What about SigmaAlpa_Coordinates part?
GetSumAcrossSigmaAlphas_withPriors <- function(LogBFs1, ModelPriors, nGammas, nSigmaAlphas) {
        WeightedSumAcrossAlphaSigmas <- matrix(0, ncol=ncol(LogBFs1), nrow=nGammas)
        for (i in 1:nGammas) {
                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
                max <- apply(LogBFs1[SigmaAlpha_Coordinates,], 2, max)
                LogBFs1[SigmaAlpha_Coordinates,] <- LogBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(LogBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(LogBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
                #20160902 CHECK_0 -- Prob: Check use of max*nSigmaAlphas at end below...point is subtracting max nSigmaAlpha number of times and then summing across those rows, so shouldn't add back that max nSigmaAlpha number of times too?
                WeightedSumAcrossAlphaSigmas[i,] <- log10(sapply(apply(ModelPriors[SigmaAlpha_Coordinates,] * apply(10^LogBFs1[SigmaAlpha_Coordinates,], c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + max*nSigmaAlphas
        }
        return(WeightedSumAcrossAlphaSigmas)
}

GetLogBFsFromData <- function(DataSources, GWASsnps, ExpectedColumnNames, SigmaAlphas, MergedDataSources, ProvidedPriors, UseFlatPrior, PruneMarginalHits, PruneMarginalHits_bpWindow, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, NminThreshold, bmassSeedValue, LogFile) {
        
	#Conducting main bmass analyses and first-level results presentation
        #~~~~~~

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Conducting main bmass analysis and first-level results formatting.", sep=""))

        if (!is.null(bmassSeedValue)) {
                set.seed(bmassSeedValue)
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- setting seed with the following value: ", bmassSeedValue, ".", sep=""))
        }
        else {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- no seed value via bmassSeedValue provided.", sep=""))
        }

        ZScoresMarginal_CommandText <- ""
        for (DataSource in DataSources) {
                if (length(strsplit(ZScoresMarginal_CommandText, "")[[1]]) == 0) {
                        ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, "cbind(MarginalHits$", DataSource, "_ZScore", sep="")
                }
                else {
                        ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, ",MarginalHits$", DataSource, "_ZScore", sep="")
                }
        }
        ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, ")", sep="")
        ZScoresMarginal <- eval(parse(text=ZScoresMarginal_CommandText))

        ZScoresMarginalNames_CommandText <- c()
        for (DataSource in DataSources) {
                ZScoresMarginalNames_CommandText <- c(ZScoresMarginalNames_CommandText, paste(DataSource, "_ZScore", sep=""))
        }
        colnames(ZScoresMarginal) <- ZScoresMarginalNames_CommandText

        NsMarginal_CommandText <- ""
        for (DataSource in DataSources) {
                if (length(strsplit(NsMarginal_CommandText, "")[[1]]) == 0) {
                        NsMarginal_CommandText <- paste(NsMarginal_CommandText, "cbind(MarginalHits$", DataSource, "_N", sep="")
                }
                else {
                        NsMarginal_CommandText <- paste(NsMarginal_CommandText, ",MarginalHits$", DataSource, "_N", sep="")
                }
        }
        NsMarginal_CommandText <- paste(NsMarginal_CommandText, ")", sep="")
        NsMarginal <- eval(parse(text=NsMarginal_CommandText))
        NsMarginal_RowMins <- apply(NsMarginal, 1, min)
        MarginalHits$Nmin <- NsMarginal_RowMins

        #20160822 20160823 CHECK_1 -- Prob: Go through use of 'do.call(rbind...etc...' and double-check logic Soln: Reminder, do.call is for applying a function to a given list of arguments. Eg the contents of do.call are treated as the full set of arguments to be used, versus say an 'apply' version where the function is applied individually to each set of arguments/vectors.
        #20160822 CHECK_0 -- Prob: Change output of Matthew's code to use logBFs vs. lbf
        #20160822 CHECK_0 -- Prob: Change output of Matthew's code to match styles developed here

        print(MarginalHits)
        print(ZScoresMarginal)

        MarginalHits_logBFs <- compute.allBFs.fromZscores(ZScoresMarginal, ZScoresCorMatrix, MarginalHits$Nmin, MarginalHits$MAF, SigmaAlphas)
        #MarginalHits_logBFs$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows
        MarginalHits_logBFs_Stacked <- do.call(rbind, MarginalHits_logBFs$lbf)

        if (!is.null(ProvidedPriors)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
                MarginalHits_logBFs$prior <- ProvidedPriors
        }

#       LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up priors.", sep=""))

#       GWASHits_EBprior <- NULL
#       FlatUnif_EBprior <- NULL
#       logBF_min
        MarginalHits_logBFs_Stacked_AvgwPrior <- NULL
        MarginalHits_logBFs_Stacked_AvgwPrior_Min <- NULL
        Priors_Used <- NULL
        if (!is.null(ProvidedPriors)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
                MarginalHits_logBFs$prior <- ProvidedPriors
                Priors_Used <- ProvidedPriors
        } else if (is.null(GWASsnps) || UseFlatPriors == TRUE) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up flat-tiered priors, GWASnps either not provided or flat prior explicitly requested.", sep=""))

                Prior_FlatUnif <- normalize(rep(c(0,MarginalHits_logBFs$prior[-1]),length(SigmaAlphas)))
                #nsigma=length(sigmaa)
                #origprior = rep(c(0,lbf$prior[-1]),nsigma)
                #origprior = normalize(origprior)

                #Prior_FlatUnif

                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and flat-tiered priors.", sep=""))

                MarginalHits_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalHits_logBFs_Stacked, Prior_FlatUnif)
                Priors_Used <- Prior_FlatUnif

                #lbf.av.origprior.glhits = lbf.av(lbf.glhits,origprior)

                #Add summary stats to marginal SNPs
                #Add SNPs x Model matrix with prior*logBFs as entries, summed (or avg'd??) across SigmaAlphas

        }
        else {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up GWAS trained priors and analyzing GWAS hits since GWASsnps provided.", sep=""))
                PreviousSNPs <- MarginalHits[MarginalHits$GWASannot==1,]
                PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalHits_logBFs_Stacked[,MarginalHits$GWASannot==1]) #Matrix of nSigmaAlphas x nSNPs

                #20160822 20160823 CHECK_1 -- Prob: Do GWAS hit analysis/work here Soln: Wrote up a few first-level GWAS hit results to get things started. Certainly will undergo further revisions down the line but fine strating point for now.

                Prior_PreviousSNPsEB <- em.priorprobs(PreviousSNPs_logBFs_Stacked, MarginalHits_logBFs$prior, 100) #Vector with nModels*nSigmaAlphas entries
                #20160823 CHECK_0: Prob -- double check use of em.priorprobs here too with runif prior starting point
                #20160823 CHECK_0: Prob -- Do multiple runs of em.priorprobs and figure out way to compare them for consistency?
                Prior_PreviousSNPsEB_check2 <- em.priorprobs(PreviousSNPs_logBFs_Stacked, MarginalHits_logBFs$prior*runif(length(MarginalHits_logBFs$prior)), 100)

                ###### do thissssss
                #
                # Run multiple EMs to check/test for convergence?
                #
                ###### do thissssss

                Prior_PreviousSNPsEB_Collapsed <- CollapseSigmaAlphasTogether(Prior_PreviousSNPsEB, length(SigmaAlphas))
                Prior_PreviousSNPsEB_check2_Collapsed <- CollapseSigmaAlphasTogether(Prior_PreviousSNPsEB_check2, length(SigmaAlphas))

                PreviousSNPs_PosteriorProbs <- posteriorprob(PreviousSNPs_logBFs_Stacked, Prior_PreviousSNPsEB) #Matrix of nModels*nSigmaAlphas x nSNPs
                PreviousSNPs_PosteriorProbs_Collapsed <- apply(PreviousSNPs_PosteriorProbs, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
#               print(PreviousSNPs_PosteriorProbs)
#               print(PreviousSNPs_PosteriorProbs_Collapsed)

###             PrevHits$BestModels <- apply(MarginalHits_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")
                PreviousSNPs$BestModel <- apply(MarginalHits_logBFs$gamma[apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, which.max),], 1, paste, collapse="_")
                PreviousSNPs$BestModel_Posterior <- apply(PreviousSNPs_PosteriorProbs_Collapsed, 2, max)
#               print(PrevHits$BestModels)
                print(PreviousSNPs)

                ## this returns a list with elements  pU pD and pI
                #marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

                Prior_PreviousSNPsEB_ModelMatrix <- cbind(MarginalHits_logBFs$gamma, Prior_PreviousSNPsEB_Collapsed)[order(Prior_PreviousSNPsEB_Collapsed, decreasing=TRUE),]
                Prior_PreviousSNPsEB_ModelMatrix <- data.frame(cbind(Prior_PreviousSNPsEB_ModelMatrix, cumsum(Prior_PreviousSNPsEB_ModelMatrix[,ncol(Prior_PreviousSNPsEB_ModelMatrix)])))
                colnames(Prior_PreviousSNPsEB_ModelMatrix) <- c(DataSources, "pValue", "Cumm_pValue")
                print(Prior_PreviousSNPsEB_ModelMatrix)

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
                PrevHits <- list()
                PrevHits$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
#               PreviousSNPs$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
#               print(Prior_PreviousSNPsEB_ModelMatrix_pValSupport)
                print(PrevHits$ModelCategories_CummpValues)

###?            MarginalHits_logBFs_Stacked_AvgwPrior <- ()
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and trained priors.", sep=""))

                MarginalHits_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalHits_logBFs_Stacked, Prior_PreviousSNPsEB)
                MarginalHits_logBFs_Stacked_AvgwPrior_Min <- min(MarginalHits_logBFs_Stacked_AvgwPrior)
                Priors_Used <- Prior_PreviousSNPsEB
        }

        #20160901 CHECK_0 -- Prob: Do something more substantive here? A better error message, or give just a warning instead? Don't exist program?
        ####Candidate For Unit Tests####
        if (is.null(MarginalHits_logBFs_Stacked_AvgwPrior)) {
                stop(Sys.time(), " -- No average log BFs were returned from method. Check if all input variables are as the method expects.")
        }

        MarginalHits$LogBFWeightedAvg <- MarginalHits_logBFs_Stacked_AvgwPrior



}





