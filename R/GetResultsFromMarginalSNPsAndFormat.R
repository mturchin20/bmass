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

#GetMeanAcrossAlphaSigmas <- function(logBFs1, nGammas, nSigmaAlphas) {
#        MeanAcrossAlphaSigmas <- matrix(0, ncol=ncol(logBFs1), nrow=nGammas)
#        for (i in 1:nGammas) {
#                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
#                max <- apply(logBFs1[SigmaAlpha_Coordinates,], 2, max)
#                logBFs1[SigmaAlpha_Coordinates,] <- logBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(logBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(logBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
#                MeanAcrossAlphaSigmas[i,] <- log10(apply(10^logBFs1[SigmaAlpha_Coordinates,], 2, mean)) + max
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
GetSumAcrossSigmaAlphas_withPriors <- function(logBFs1, ModelPriors, nGammas, nSigmaAlphas) {
        WeightedSumAcrossAlphaSigmas <- matrix(0, ncol=ncol(logBFs1), nrow=nGammas)
        for (i in 1:nGammas) {
                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
                max <- apply(logBFs1[SigmaAlpha_Coordinates,], 2, max)
                logBFs1[SigmaAlpha_Coordinates,] <- logBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(logBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(logBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
                #20160902 CHECK_0 -- Prob: Check use of max*nSigmaAlphas at end below...point is subtracting max nSigmaAlpha number of times and then summing across those rows, so shouldn't add back that max nSigmaAlpha number of times too?
                WeightedSumAcrossAlphaSigmas[i,] <- log10(sapply(apply(ModelPriors[SigmaAlpha_Coordinates,] * apply(10^logBFs1[SigmaAlpha_Coordinates,], c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + max*nSigmaAlphas
        }
        return(WeightedSumAcrossAlphaSigmas)
}

GetLogBFsFromData <- function(DataSources, MarginalSNPs, ZScoresCorMatrix, SigmaAlphas, LogFile) {
        
	#Conducting main bmass analyses and first-level results presentation
        #~~~~~~

	PreviousSNPs <- list()

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Conducting main bmass analysis and first-level results formatting.", sep=""))

        ZScoresMarginal_CommandText <- ""
        for (DataSource in DataSources) {
                if (length(strsplit(ZScoresMarginal_CommandText, "")[[1]]) == 0) {
                        ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, "cbind(MarginalSNPs$SNPs$", DataSource, "_ZScore", sep="")
                }
                else {
                        ZScoresMarginal_CommandText <- paste(ZScoresMarginal_CommandText, ",MarginalSNPs$SNPs$", DataSource, "_ZScore", sep="")
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
                        NsMarginal_CommandText <- paste(NsMarginal_CommandText, "cbind(MarginalSNPs$SNPs$", DataSource, "_N", sep="")
                }
                else {
                        NsMarginal_CommandText <- paste(NsMarginal_CommandText, ",MarginalSNPs$SNPs$", DataSource, "_N", sep="")
                }
        }
        NsMarginal_CommandText <- paste(NsMarginal_CommandText, ")", sep="")
        NsMarginal <- eval(parse(text=NsMarginal_CommandText))
        NsMarginal_RowMins <- apply(NsMarginal, 1, min)
        MarginalSNPs$SNPs$Nmin <- NsMarginal_RowMins

        #20160822 20160823 CHECK_1 -- Prob: Go through use of 'do.call(rbind...etc...' and double-check logic Soln: Reminder, do.call is for applying a function to a given list of arguments. Eg the contents of do.call are treated as the full set of arguments to be used, versus say an 'apply' version where the function is applied individually to each set of arguments/vectors.
        #20160822 CHECK_0 -- Prob: Change output of Matthew's code to use logBFs vs. lbf
        #20160822 CHECK_0 -- Prob: Change output of Matthew's code to match styles developed here

#       print(MarginalSNPs)
#       print(ZScoresMarginal)

        MarginalSNPs_logBFs <- compute.allBFs.fromZscores(ZScoresMarginal, ZScoresCorMatrix, MarginalSNPs$SNPs$Nmin, MarginalSNPs$SNPs$MAF, SigmaAlphas)
        MarginalSNPs_logBFs_Stacked <- do.call(rbind, MarginalSNPs_logBFs$lbf)
	MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked

	return(list(MarginalSNPs=MarginalSNPs, Models=MarginalSNPs_logBFs$gamma, ModelPriors=MarginalSNPs_logBFs$prior, LogFile=LogFile))

}

DetermineAndApplyPriors <- function(DataSources, MarginalSNPs, GWASsnps, SigmaAlphas, Models, ModelPriors, ProvidedPriors, UseFlatPriors, bmassSeedValue, LogFile) {

	MarginalSNPs_logBFs_Stacked <- MarginalSNPs$logBFs	
	MarginalSNPs_logBFs_Stacked_AvgwPrior <- NULL
        ModelPriors_Used <- ModelPriors
       
	PreviousSNPs <- list()
        PreviousSNPs_logBFs_Stacked_AvgwPrior_Min <- NULL
 
	if (!is.null(ProvidedPriors)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, ProvidedPriors)
                ModelPriors_Used <- ProvidedPriors
        } else if (is.null(GWASsnps) || UseFlatPriors == TRUE) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up flat-tiered priors, GWASnps either not provided or flat prior explicitly requested.", sep=""))

                #20160930 CHECK_0 -- Prob: Go over this use again of 0 to start the original prior setup w/ Matthew
		Prior_FlatUnif <- normalize(rep(c(0,ModelPriors[-1]),length(SigmaAlphas)))
                #nsigma=length(sigmaa)
                #origprior = rep(c(0,lbf$prior[-1]),nsigma)
                #origprior = normalize(origprior)

                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and flat-tiered priors.", sep=""))

                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, Prior_FlatUnif)
                ModelPriors_Used <- Prior_FlatUnif

		#20160930 CHECK_0 -- Prob: Do this below? Had no CHECK_0 note until this one so maybe just passing reminder/note that was or isn't meant to be taken care of eventually?
                #Add summary stats to marginal SNPs
                #Add SNPs x Model matrix with prior*logBFs as entries, summed (or avg'd??) across SigmaAlphas
        }
        else {
	
		if (!is.null(bmassSeedValue)) {
			set.seed(bmassSeedValue)
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- setting seed with the following value: ", bmassSeedValue, ".", sep=""))
		}
		else {
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- no seed value via bmassSeedValue provided.", sep=""))
		}
		
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up GWAS trained priors and analyzing GWAS hits since GWASsnps provided.", sep=""))
                
		PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1]) #Matrix of nSigmaAlphas x nSNPs

                #20160822 20160823 CHECK_1 -- Prob: Do GWAS hit analysis/work here Soln: Wrote up a few first-level GWAS hit results to get things started. Certainly will undergo further revisions down the line but fine strating point for now.

                Prior_PreviousSNPsEB <- em.priorprobs(PreviousSNPs_logBFs_Stacked, ModelPriors, 100) #Vector with nModels*nSigmaAlphas entries
                #20160823 CHECK_0: Prob -- double check use of em.priorprobs here too with runif prior starting point
                #20160823 CHECK_0: Prob -- Do multiple runs of em.priorprobs and figure out way to compare them for consistency?
                Prior_PreviousSNPsEB_check2 <- em.priorprobs(PreviousSNPs_logBFs_Stacked, ModelPriors*runif(length(ModelPriors)), 100)

                ###### do thissssss
                #
                # Run multiple EMs to check/test for convergence?
                #
                ###### do thissssss

                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, Prior_PreviousSNPsEB)
                ModelPriors_Used <- Prior_PreviousSNPsEB
                
		PreviousSNPs$logBFs <- PreviousSNPs_logBFs_Stacked
        }

        #20160901 CHECK_0 -- Prob: Do something more substantive here? A better error message, or give just a warning instead? Don't exist program?
        ####Candidate For Unit Tests####
        if (is.null(MarginalSNPs_logBFs_Stacked_AvgwPrior)) {
                stop(Sys.time(), " -- No average log BFs were returned from method. Check if all input variables are as the method expects.")
        }
        
	MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked
        MarginalSNPs$SNPs$logBFWeightedAvg <- MarginalSNPs_logBFs_Stacked_AvgwPrior
		
	PreviousSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1,]
        if (dim(PreviousSNPs$SNPs)[1] > 0) {
		PreviousSNPs_logBFs_Stacked_AvgwPrior_Min <- min(PreviousSNPs$SNPs$logBFWeightedAvg)
	}

	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, ModelPriors=ModelPriors_Used, GWASlogBFMinThreshold=PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, LogFile=LogFile))

}


FinalizeAndFormatResults <- function(DataSources, MarginalSNPs, PreviousSNPs, GWASsnps, PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, SigmaAlphas, Models, ModelPriors, NminThreshold, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, LogFile) {

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and trained priors.", sep=""))
	
        NewSNPs <- list()
	MarginalSNPs_logBFs_Stacked <- MarginalSNPs$logBF
	PreviousSNPs_logBFs_Stacked <- PreviousSNPs$logBF

        #Pruning marginal hits by logBFWeightedAvg if requested
        if (PruneMarginalSNPs == TRUE) {
                #20160901 CHECK_0 -- Prob: Go over indepthits function, rewrite, or just lightly edit? redo names, double-check functionality? def get some unit testing in there
                MarginalSNPs_PrunedList <- indephits(MarginalSNPs$SNPs$logBFWeightedAvg, MarginalSNPs$SNPs$Chr, MarginalSNPs$SNPs$BP, T=PruneMarginalSNPs_bpWindow)
                MarginalSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs_PrunedList==1,]
                MarginalSNPs_logBFs_Stacked <- MarginalSNPs_logBFs_Stacked[,MarginalSNPs_PrunedList==1]
        }

        #Summing models over all values of SigmaAlphas, weighted by ModelPriors
	MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed <- GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(ModelPriors, ncol(MarginalSNPs_logBFs_Stacked)), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), nrow(Models), length(SigmaAlphas))
        MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed <- cbind(Models, MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed)
        colnames(MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed) <- c(DataSources, MarginalSNPs$SNPs$ChrBP)
	MarginalSNPs$logBF <- MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed
        
	PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed <- GetSumAcrossSigmaAlphas_withPriors(PreviousSNPs_logBFs_Stacked, matrix(rep(ModelPriors, ncol(PreviousSNPs_logBFs_Stacked)), ncol=ncol(PreviousSNPs_logBFs_Stacked), byrow=FALSE), nrow(Models), length(SigmaAlphas))
        PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed <- cbind(Models, PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed)
        colnames(PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed) <- c(DataSources, PreviousSNPs$SNPs$ChrBP)
	PreviousSNPs$logBF <- PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed

        #Preparing posterior probabilities, Gammas x SNPs
        MarginalSNPs_logBFs_Stacked_Posteriors <- posteriorprob(MarginalSNPs_logBFs_Stacked, ModelPriors)
        MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed <- apply(MarginalSNPs_logBFs_Stacked_Posteriors, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas))
        MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed <- cbind(Models, MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed)
        colnames(MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed) <- c(DataSources, MarginalSNPs$SNPs$ChrBP)
	MarginalSNPs$Posteriors <- MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed

	PreviousSNPs_logBFs_Stacked_Posteriors <- posteriorprob(PreviousSNPs_logBFs_Stacked, ModelPriors) #Matrix of nModels*nSigmaAlphas x nSNPs
        PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed <- apply(PreviousSNPs_logBFs_Stacked_Posteriors, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
        PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed <- cbind(Models, PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed)
        colnames(PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed) <- c(DataSources, PreviousSNPs$SNPs$ChrBP)
	PreviousSNPs$Posteriors <- PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed

        #20160905 20161005 CHECK_1 -- Prob: Convert either 'MarginalSNPs_logBFs_Stacked_PosteriorProbabilities' to 'PosteriorProbs' or 'PreviousSNPs_PosteriorProbs' to 'PreviousSNPs_PosteriorProbabilities' Soln: Actually deicded to change both terms to 'Posteriors', so neither 'PosteriorProbabilities' or 'PosteriorProbs' was used. 
        #20160905 201609** CHECK_1 -- Prob: Change 'MarginalSNPs' to 'MarginalSNPs' probably? Soln: I made this change as evidenced by the aforementioned names are the same -- did a file-wide substitution of 'MarginalHits' to 'MarginalSNPs' so that's probably what was being referenced here. Didn't make the note here when I made the change so I don't recall the exact date I did this, but it was sometime before 20161005 and likely mid/late September.

        #20160905 20161005 CHECK_1 -- Prob: Move this (NewSNPs <- NULL) below intialization section to proper beginning of .R file code as necessary once change/reorganization occurs Soln: Moved it to top of this function block

        #Determining new hits if GWASsnps were provided to determine minimum MarginalSNPs_logBFs_Stacked_AvgwPrior value threshold
        #20160902 CHECK_0 -- Prob: Check that PreviousSNPs_logBFs_Stacked_AvgwPrior_Min is non null here?
        if (!is.null(GWASsnps)) {
                if (is.null(PreviousSNPs_logBFs_Stacked_AvgwPrior_Min)) {
                        stop(Sys.time(), " -- PreviousSNPs_logBFs_Stacked_AvgwPrior_Min is NULL despite GWASsnps being provided. Unexpected error.")
                }
                NewSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot == 0 & MarginalSNPs$SNPs$logBFWeightedAvg >= PreviousSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$SNPs$Nmin >= NminThreshold,]
                NewSNPs$logBF <- MarginalSNPs$logBF[,MarginalSNPs$SNPs$GWASannot == 0 & MarginalSNPs$SNPs$logBFWeightedAvg >= PreviousSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$SNPs$Nmin >= NminThreshold]
                NewSNPs$Posteriors <- MarginalSNPs$Posteriors[,MarginalSNPs$SNPs$GWASannot == 0 & MarginalSNPs$SNPs$logBFWeightedAvg >= PreviousSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$SNPs$Nmin >= NminThreshold]
        }

        ####Candidate For Unit Tests####
        #print(dim(matrix(rep(ModelPriors, ncol(MarginalSNPs_logBFs_Stacked)), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE)))
        
	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, NewSNPs=NewSNPs, LogFile=LogFile))

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
#	PrevHits <- list()
#	PrevHits$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
##	PreviousSNPs$ModelCategories_CummpValues <- Prior_PreviousSNPsEB_ModelMatrix_pValSupport
##	print(PrevHits$ModelCategories_CummpValues)



	
	
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


