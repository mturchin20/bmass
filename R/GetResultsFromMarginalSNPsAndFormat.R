#This function expects columns from nSigmaAlphas stacked Model x SNP matrices of posterior probabilities, such that a single column (ie single SNP) is converted to a Model x nSigmaAlphas matrix and summed across rows for a single posterior probability per Model (eg the 'marginal' of across all sigma_alphas)  
CollapseSigmaAlphasTogether <- function (inputValues1, nSigmaAlphas) {
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

####Candidates For Unit Tests####
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

#This function is expecting nSigmaAlphas Model x SNP matrices of logBFs stacked ontop of one another. ModelPriors_Matrix is a matrix containing the vector of model priors (ModelPriors = one set of model priors * nSigmaAlphas) replicated as a column for each snp. 
GetSumAcrossSigmaAlphas_withPriors <- function(logBFs1, ModelPriors_Matrix, nGammas, nSigmaAlphas) {
	WeightedSumAcrossAlphaSigmas <- matrix(0, ncol=ncol(logBFs1), nrow=nGammas)
        for (i in 1:nGammas) {
                SigmaAlpha_Coordinates <- seq.int(from=i, by=nGammas, length.out=nSigmaAlphas)
                max <- apply(logBFs1[SigmaAlpha_Coordinates,], 2, max)
                logBFs1[SigmaAlpha_Coordinates,] <- logBFs1[SigmaAlpha_Coordinates,] - matrix(max, nrow=nrow(logBFs1[SigmaAlpha_Coordinates,]), ncol=ncol(logBFs1[SigmaAlpha_Coordinates,]), byrow=TRUE)
		WeightedSumAcrossAlphaSigmas[i,] <- log10(sapply(apply(ModelPriors_Matrix[SigmaAlpha_Coordinates,] * apply(10^logBFs1[SigmaAlpha_Coordinates,], c(1,2), CheckForAndReplaceOnes), 2, sum), CheckForAndReplaceZeroes)) + max
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

        MarginalSNPs_logBFs <- compute.allBFs.fromZscores(ZScoresMarginal, ZScoresCorMatrix, MarginalSNPs$SNPs$Nmin, MarginalSNPs$SNPs$MAF, SigmaAlphas)
        MarginalSNPs_logBFs_Stacked <- do.call(rbind, MarginalSNPs_logBFs$lbf)
	MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked

	return(list(MarginalSNPs=MarginalSNPs, Models=MarginalSNPs_logBFs$gamma, ModelPriors=MarginalSNPs_logBFs$prior, LogFile=LogFile))

}

DetermineAndApplyPriors <- function(DataSources, MarginalSNPs, GWASsnps, SigmaAlphas, Models, ModelPriors, ProvidedPriors, UseFlatPriors, GWASThreshFlag, GWASThreshValue, bmassSeedValue, LogFile) {

	MarginalSNPs_logBFs_Stacked <- MarginalSNPs$logBFs	
	MarginalSNPs_logBFs_Stacked_AvgwPrior <- NULL
        ModelPriors_Used <- ModelPriors
		
	PreviousSNPs <- list()
        PreviousSNPs_logBFs_Stacked_AvgwPrior_Min <- NULL
	ZScoreHitFlag1 <- c()
	if (GWASThreshFlag) {
		ZScoreHitFlag1 <- rep(0, nrow(MarginalSNPs$SNPs))
		ZScoreHitFlag1[2*pnorm(apply(abs(MarginalSNPs$SNPs[,grep("ZScore", colnames(MarginalSNPs$SNPs))]),1,max),0,1,lower.tail=FALSE) < GWASThreshValue] <- 1
	}

	if (!is.null(ProvidedPriors)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- ProvidedPriors is not NULL, replacing original priors with submitted values.", sep=""))
                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, ProvidedPriors)
                ModelPriors_Used <- ProvidedPriors
		PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1]) #Matrix of nSigmaAlphas x nSNPs
		PreviousSNPs$logBFs <- PreviousSNPs_logBFs_Stacked
        } else if (is.null(GWASsnps) && UseFlatPriors == TRUE) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up flat-tiered priors, GWASnps either not provided or flat prior explicitly requested.", sep=""))

		Prior_FlatUnif <- normalize(rep(c(0,ModelPriors[-1]),length(SigmaAlphas)))

                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and flat-tiered priors.", sep=""))

                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, Prior_FlatUnif)
                ModelPriors_Used <- Prior_FlatUnif

        } else if (is.null(GWASsnps) && UseFlatPriors == FALSE) {
		#20171017 CHECK_0 -- Prob: Throw an error/fail thing here
		PH <- 1
	} else {
	
		if (!is.null(bmassSeedValue)) {
			set.seed(bmassSeedValue)
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- setting seed with the following value: ", bmassSeedValue, ".", sep=""))
		}
		else {
			LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- no seed value via bmassSeedValue provided.", sep=""))
		}
		
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Setting up GWAS trained priors and analyzing GWAS hits since GWASsnps provided.", sep=""))
               
		PreviousSNPs_logBFs_Stacked <- c() #Will be matrix of nSigmaAlphas x nSNPs
		if (GWASThreshFlag) {
			PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1 & ZScoreHitFlag1==1]) 
		} 
		else if (!GWASThreshFlag) {
			PreviousSNPs_logBFs_Stacked <- as.matrix(MarginalSNPs_logBFs_Stacked[,MarginalSNPs$SNPs$GWASannot==1]) 
		} 
		else {
			#20171017 CHECK_0 -- Prob: Throw an error/fail thing here
			PH <- 1
		}

                Prior_PreviousSNPsEB <- em.priorprobs(PreviousSNPs_logBFs_Stacked, ModelPriors, 100) #Vector with nModels*nSigmaAlphas entries
                Prior_PreviousSNPsEB_check2 <- em.priorprobs(PreviousSNPs_logBFs_Stacked, ModelPriors*runif(length(ModelPriors)), 100)

                MarginalSNPs_logBFs_Stacked_AvgwPrior <- lbf.av(MarginalSNPs_logBFs_Stacked, Prior_PreviousSNPsEB)
                ModelPriors_Used <- Prior_PreviousSNPsEB
                
		PreviousSNPs$logBFs <- PreviousSNPs_logBFs_Stacked
        }

        if (is.null(MarginalSNPs_logBFs_Stacked_AvgwPrior)) {
                stop(Sys.time(), " -- No average log BFs were returned from method. Check if all input variables are as the method expects.")
        }
        
	MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked
        MarginalSNPs$SNPs$logBFWeightedAvg <- MarginalSNPs_logBFs_Stacked_AvgwPrior
		
	if (!is.null(GWASsnps)) {
		if (GWASThreshFlag) {
			PreviousSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1 & ZScoreHitFlag1==1,]
			PreviousSNPs$DontPassSNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1 & ZScoreHitFlag1==0,]
		} else if (!GWASThreshFlag) {
			PreviousSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot==1,]
		} else {
			#20171017 CHECK_0 -- Prob: Throw an error/fail thing here
			PH <- 1
		}	
	
		if (dim(PreviousSNPs$SNPs)[1] > 0) {
			PreviousSNPs_logBFs_Stacked_AvgwPrior_Min <- min(PreviousSNPs$SNPs$logBFWeightedAvg)
		}
	}

	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, ModelPriors=ModelPriors_Used, GWASlogBFMinThreshold=PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, LogFile=LogFile))

}


FinalizeAndFormatResults <- function(DataSources, MarginalSNPs, PreviousSNPs, GWASsnps, PreviousSNPs_logBFs_Stacked_AvgwPrior_Min, SigmaAlphas, Models, ModelPriors, NminThreshold, PruneMarginalSNPs, PruneMarginalSNPs_bpWindow, LogFile) {

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Identifying potential new hits based on average log BFs and trained priors.", sep=""))
	
        NewSNPs <- list()
	MarginalSNPs_logBFs_Stacked <- MarginalSNPs$logBFs

        if (PruneMarginalSNPs == TRUE) {
                #20160901 20180829 CHECK_1 -- Prob: Go over indepthits function, rewrite, or just lightly edit? redo names, double-check functionality? def get some unit testing in there Soln: not going to rewrite (for now), but pushing this into the rest of the 'to create unit tests for' notes I'm taking at the moment
                MarginalSNPs_PrunedList <- indephits(MarginalSNPs$SNPs$logBFWeightedAvg, MarginalSNPs$SNPs$Chr, MarginalSNPs$SNPs$BP, T=PruneMarginalSNPs_bpWindow)
                MarginalSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs_PrunedList==1,]
                MarginalSNPs_logBFs_Stacked <- MarginalSNPs_logBFs_Stacked[,MarginalSNPs_PrunedList==1]
        }

	#Summing models over all values of SigmaAlphas, weighted by ModelPriors
	MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed <- GetSumAcrossSigmaAlphas_withPriors(MarginalSNPs_logBFs_Stacked, matrix(rep(ModelPriors, ncol(MarginalSNPs_logBFs_Stacked)), nrow=length(ModelPriors), ncol=ncol(MarginalSNPs_logBFs_Stacked), byrow=FALSE), nrow(Models), length(SigmaAlphas))
	MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed <- cbind(Models, MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed)
	colnames(MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed) <- c(DataSources, MarginalSNPs$SNPs$ChrBP)
	MarginalSNPs$logBFs <- MarginalSNPs_logBFs_Stacked_SigmaAlphasSummed

        #Preparing posterior probabilities, Gammas x SNPs
        MarginalSNPs_logBFs_Stacked_Posteriors <- posteriorprob(MarginalSNPs_logBFs_Stacked, ModelPriors)
        MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed <- apply(MarginalSNPs_logBFs_Stacked_Posteriors, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas))
        MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed <- cbind(Models, MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed)
        colnames(MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed) <- c(DataSources, MarginalSNPs$SNPs$ChrBP)
	MarginalSNPs$Posteriors <- MarginalSNPs_logBFs_Stacked_Posteriors_Collapsed

	if (!is.null(GWASsnps)) {
		PreviousSNPs_logBFs_Stacked <- PreviousSNPs$logBFs

		PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed <- GetSumAcrossSigmaAlphas_withPriors(PreviousSNPs_logBFs_Stacked, matrix(rep(ModelPriors, ncol(PreviousSNPs_logBFs_Stacked)), nrow=length(ModelPriors), ncol=ncol(PreviousSNPs_logBFs_Stacked), byrow=FALSE), nrow(Models), length(SigmaAlphas))
		PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed <- cbind(Models, PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed)
		colnames(PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed) <- c(DataSources, PreviousSNPs$SNPs$ChrBP)
		PreviousSNPs$logBFs <- PreviousSNPs_logBFs_Stacked_SigmaAlphasSummed
	
		PreviousSNPs_logBFs_Stacked_Posteriors <- posteriorprob(PreviousSNPs_logBFs_Stacked, ModelPriors) #Matrix of nModels*nSigmaAlphas x nSNPs
        	PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed <- apply(PreviousSNPs_logBFs_Stacked_Posteriors, 2, CollapseSigmaAlphasTogether, nSigmaAlphas=length(SigmaAlphas)) #Matrix of nModels x nSNPs
        	PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed <- cbind(Models, PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed)
        	colnames(PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed) <- c(DataSources, PreviousSNPs$SNPs$ChrBP)
		PreviousSNPs$Posteriors <- PreviousSNPs_logBFs_Stacked_Posteriors_Collapsed
	}

        #Determining new hits if GWASsnps were provided to determine minimum MarginalSNPs_logBFs_Stacked_AvgwPrior value threshold
        if (!is.null(GWASsnps)) {
                if (is.null(PreviousSNPs_logBFs_Stacked_AvgwPrior_Min)) {
                        stop(Sys.time(), " -- PreviousSNPs_logBFs_Stacked_AvgwPrior_Min is NULL despite GWASsnps being provided. Unexpected error.")
                }
                NewSNPs$SNPs <- MarginalSNPs$SNPs[MarginalSNPs$SNPs$GWASannot == 0 & MarginalSNPs$SNPs$logBFWeightedAvg >= PreviousSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$SNPs$Nmin >= NminThreshold,]
		NewSNPs$logBFs <- cbind(MarginalSNPs$logBFs[,1:length(DataSources)], MarginalSNPs$logBFs[,(length(DataSources)+1):ncol(MarginalSNPs$logBFs)][,MarginalSNPs$SNPs$GWASannot == 0 & MarginalSNPs$SNPs$logBFWeightedAvg >= PreviousSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$SNPs$Nmin >= NminThreshold])
                NewSNPs$Posteriors <- cbind(MarginalSNPs$Posteriors[,1:length(DataSources)], MarginalSNPs$Posteriors[,(length(DataSources)+1):ncol(MarginalSNPs$Posteriors)][,MarginalSNPs$SNPs$GWASannot == 0 & MarginalSNPs$SNPs$logBFWeightedAvg >= PreviousSNPs_logBFs_Stacked_AvgwPrior_Min & MarginalSNPs$SNPs$Nmin >= NminThreshold])
        }

	return(list(MarginalSNPs=MarginalSNPs, PreviousSNPs=PreviousSNPs, NewSNPs=NewSNPs, LogFile=LogFile))
}
