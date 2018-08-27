
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

CountModelClasses <- function(ModelEntries) { 
	Count0Class <- 0 
	Count1Class <- 0 
	Count2Class <- 0 
	for (j in 1:length(ModelEntries)) { 
		if (ModelEntries[j] == 0) { 
			Count0Class <- Count0Class + 1
		} else if (ModelEntries[j] == 1) { 
			Count1Class <- Count1Class + 1 
		} else if (ModelEntries[j] == 2) { 
			Count2Class <- Count2Class + 1 
		} else {
			PH1 <- 1
		}
	} 
	return(c(Count0Class, Count1Class, Count2Class))
}

#20170219 NOTE -- Unsure if going to really be using the below function eventually. Might hack together components separately and then reformulate a different 'main' function that bring them together, all outside of the below function here, than resue or continuing to operate on this one. Feel like this was created as a placeholder and not touched for a while, so contains elements that don't follow where/how the package has gone elsewhere
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

GetModelPriorMatrix <- function (DataSources, Models, ModelPriors, SigmaAlphas, LogFile) {

	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff1.", sep=""))

	ModelPriors.Collapsed <- CollapseSigmaAlphasTogether(ModelPriors, length(SigmaAlphas))

	ModelPriorMatrix <- cbind(Models,ModelPriors.Collapsed)[order(ModelPriors.Collapsed, decreasing=TRUE),]
	ModelPriorMatrix <- data.frame(cbind(ModelPriorMatrix, cumsum(ModelPriorMatrix[,(length(DataSources)+1)])), order(ModelPriors.Collapsed, decreasing=TRUE))
	#20170219 20170219 CHECK_1 -- Prob: Preserve order another way? Soln: Decided this is fine because if move to do a sort of 'sapply(which...' setup also need to import `Models` since neither `ListSNPs` or `ModelPriorMatrix` below has that original information otherwise
	colnames(ModelPriorMatrix) <- c(DataSources, "pValue", "Cumm_pValue", "OrigOrder")

	return(list(ModelPriorMatrix=ModelPriorMatrix, LogFile=LogFile))

}

#20170220 NOTE -- ModelMatrix currently being treated as independent few lines of code outside of bmass package, so not doing it inside here yet...
GetTopModelsPerSNPViaPosteriors <- function (DataSources, ListSNPs, ModelPriorMatrix, LogFile) {

	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff2.", sep=""))
	
	ListSNPs.Posteriors <- cbind(ListSNPs$Posteriors[,(length(DataSources)+1):ncol(ListSNPs$Posteriors)])[ModelPriorMatrix[,ncol(ModelPriorMatrix)],]
	
#	ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]

	ModelPriorMatrix.Paste <- apply(ModelPriorMatrix[,1:length(DataSources)], 1, function(x) paste(x, collapse="_"))
	PerSNPTopModels <- cbind(ModelPriorMatrix.Paste[apply(ListSNPs.Posteriors, 2, which.max)], apply(ListSNPs.Posteriors, 2, max))
	SummaryOfTopModels <- c()
	
	for (i in unique(PerSNPTopModels[,1])) {
	        PerSNPTopModels.Subset <- matrix(PerSNPTopModels[PerSNPTopModels[,1]==i,], nrow=length(which(PerSNPTopModels[,1]==i)))
	        SummaryOfTopModels <- rbind(SummaryOfTopModels, cbind(i, nrow(PerSNPTopModels.Subset), round(mean(as.numeric(PerSNPTopModels.Subset[,2])), digits=3), ModelPriorMatrix[ModelPriorMatrix.Paste == i, length(DataSources)+1]))
	}

	if (nrow(SummaryOfTopModels) > 1) {
		SummaryOfTopModels <- SummaryOfTopModels[order(as.numeric(SummaryOfTopModels[,2]), decreasing=TRUE),]
	}
	colnames(SummaryOfTopModels) <- c(paste(DataSources, collapse="_"), "n", "MeanPosterior", "OriginalPrior")
	ListSNPs$TopModels <- SummaryOfTopModels

        return(list(ListSNPs=ListSNPs, LogFile=LogFile))

}


GetMarginalPosteriors <- function (DataSources, ListSNPs, Models, SigmaAlphas, LogFile) {
	
	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff3.", sep=""))
        
#	ListSNPs.MarginalPosteriors <- marginal.postprobs(ListSNPs$Posteriors[,(length(DataSources)+1):ncol(ListSNPs$Posteriors)], Models, length(SigmaAlphas))	
	ListSNPs.MarginalPosteriors <- marginal.postprobs(ListSNPs$Posteriors[,(length(DataSources)+1):ncol(ListSNPs$Posteriors)], Models, 1)	
	rownames(ListSNPs.MarginalPosteriors$pU) <- DataSources 

	ListSNPs$Marginals <- ListSNPs.MarginalPosteriors

	return(list(ListSNPs=ListSNPs, LogFile=LogFile))

}

#20170221 CHECK_0 -- Prob: Make separate version with direction?
#20170222 CHECK_0 -- Prob: keep as one function?
#20170222 CHECK_0 -- Prob: Calling ggplot on own doesn't work (though plot() does), so hack of 'print(ggplot...' works but why? Fix this, right?
#20170221 CHECK_0 -- Prob: Control/offer choice on printing out reordered marginal results?
#plot options?
#PlotMarginalPosteriors <- function (DataSources, ListSNPs, Marginal, PrintMarginals=1) {
PlotMarginalPosteriors <- function (DataSources, ListSNPs, Marginal, PrintMarginals=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="white", highcolor="steelblue") {

#	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff4.", sep=""))

	Marginal.Position <- -9

	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}

	ListSNPs.MarginalPosteriors <- t(ListSNPs$Marginals[[Marginal.Position]])
	ListSNPs.MarginalPosteriors <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.Melted <- melt(ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.Melted$Marker, levels=ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors <- ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,]

#	print(ggplot(ListSNPs.MarginalPosteriors.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = lowcolor, high=highcolor) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
#20170622 NOTE -- included the `limits` addition and the movement to a returning of the plot in order to help with paper/manuscript/chapter stuff -- not meant to be included for final package. `limits` may be worthwhile to keep as an option though.
	print(ggplot(ListSNPs.MarginalPosteriors.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(limits = c(-0,1), low = lowcolor, high=highcolor) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))

	if (PrintMarginals == 0) {
		return (PH <- 1)
	} else if (PrintMarginals == 1) {
		return(ListSNPs.MarginalPosteriors)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}


}

#20170222 NOTE -- Below has not been tested yet
PlotMarginalPosteriors.withDirection <- function (DataSources, ListSNPs, Marginal, PrintMarginals=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="darkred", midcolor="white", highcolor="steelblue") {

#	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff4.", sep=""))

	Marginal.Position <- -9

	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}

	ListSNPs.MarginalPosteriors <- t(ListSNPs$Marginals[[Marginal.Position]])
	ListSNPs.MarginalPosteriors <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors) <- c("Marker", DataSources)	

	Z <- c()
	for (i in DataSources) {
		eval(parse(text=paste("Z <- cbind(Z, ListSNPs$SNPs$", i, "_ZScore)", sep="")))
	}

	Z.Directions <- apply(Z, c(1,2), function(x) { y <- 0; if (x > 0) { y <- 1; } else if (x < 0) { y <- -1; } else if (x == 0) { y <- 0; } else { y <- NA; }; return(y);})

	ListSNPs.MarginalPosteriors.wDirection <- ListSNPs.MarginalPosteriors
	ListSNPs.MarginalPosteriors.wDirection[,2:(length(DataSources)+1)] <- ListSNPs.MarginalPosteriors.wDirection[,2:(length(DataSources)+1)] * Z.Directions

	ListSNPs.MarginalPosteriors.Melted <- melt(ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.Melted$Marker, levels=ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.wDirection.Melted <- melt(ListSNPs.MarginalPosteriors.wDirection[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.wDirection.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.wDirection.Melted$Marker, levels=ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors <- ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,]

	print(ggplot(ListSNPs.MarginalPosteriors.wDirection.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))

	if (PrintMarginals == 0) {
		return (PH <- 1)
	} else if (PrintMarginals == 1) {
		return(ListSNPs.MarginalPosteriors)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}

}

#20170301 NOTE -- Below has not been tested yet
PlotZScores.OrderedByClusteredMarginals <- function (DataSources, ListSNPs, Marginal, PrintZScores=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="darkred", midcolor="white", highcolor="steelblue") {

	Marginal.Position <- -9

	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}

	ListSNPs.MarginalPosteriors <- t(ListSNPs$Marginals[[Marginal.Position]])
	ListSNPs.MarginalPosteriors <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors) <- c("Marker", DataSources)	

	ZScores <- data.frame(ListSNPs$SNPs$Marker)
	for (i in DataSources) {
		eval(parse(text=paste("ZScores <- cbind(ZScores, ListSNPs$SNPs$", i, "_ZScore)", sep="")))
	}
	colnames(ZScores) <- c("Marker", DataSources)	

	ZScores.Melted <- melt(ZScores[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,])
	ZScores.Melted$Marker <- factor(ZScores.Melted$Marker, levels=ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,1])

	ggplot(ZScores.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main)

	print(ggplot(ZScores.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))

	ZScores <- ZScores[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,]
	
	if (PrintZScores == 0) {
		return (PH <- 1)
	} else if (PrintZScores == 1) {
		return(ZScores)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}

}

#####
#
# 20171019 NOTE -- Below code/functions meant for manuscript production, not final release of the R package (possibly make this a separate .R file?)
#
#####


#20170622 NOTE -- added/made this function below
PlotMarginalPosteriors.ForManuscript <- function (DataSources, ListSNPs, Marginal, PrintMarginals=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="white", highcolor="steelblue") {

#	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff4.", sep=""))

	Marginal.Position <- -9

	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}

	ListSNPs.MarginalPosteriors.D <- t(ListSNPs$Marginals[[2]])
	ListSNPs.MarginalPosteriors.D <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors.D))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors.D[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors.D[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors.D) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.D.Melted <- melt(ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.D.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.D.Melted$Marker, levels=ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.I <- t(ListSNPs$Marginals[[3]])
	ListSNPs.MarginalPosteriors.I <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors.I))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors.I[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors.I[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors.I) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.I.Melted <- melt(ListSNPs.MarginalPosteriors.I[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.I.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.I.Melted$Marker, levels=ListSNPs.MarginalPosteriors.I[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.U <- t(ListSNPs$Marginals[[1]])
	ListSNPs.MarginalPosteriors.U <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors.U))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors.U[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors.U[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors.U) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.U.Melted <- melt(ListSNPs.MarginalPosteriors.U[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.U.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.U.Melted$Marker, levels=ListSNPs.MarginalPosteriors.U[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

#20170622 NOTE -- included the `limits` addition and the movement to a returning of the plot in order to help with paper/manuscript/chapter stuff -- not meant to be included for final package. `limits` may be worthwhile to keep as an option though.
#	print(ggplot(ListSNPs.MarginalPosteriors.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(limits = c(-0,1), low = lowcolor, high=highcolor) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))

	if (PrintMarginals == 0) {
		if (Marginal == "U") {
			return(ggplot(ListSNPs.MarginalPosteriors.U.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(limits = c(-0,1), low = lowcolor, high=highcolor) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
		} else if (Marginal == "D") {
			return(ggplot(ListSNPs.MarginalPosteriors.D.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(limits = c(-0,1), low = lowcolor, high=highcolor) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
		} else if (Marginal == "I") {
			return(ggplot(ListSNPs.MarginalPosteriors.I.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(limits = c(-0,1), low = lowcolor, high=highcolor) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
		} else {
			return(PH <- 1)
		}
	} else if (PrintMarginals == 1) {
		return(ListSNPs.MarginalPosteriors)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}

}


#20170622 NOTE -- added/made this function below
PlotMarginalPosteriors.withDirection.ForManuscript <- function (DataSources, ListSNPs, Marginal, PrintMarginals=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="darkred", midcolor="white", highcolor="steelblue") {

#	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff4.", sep=""))

	Marginal.Position <- -9

	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}

	ListSNPs.MarginalPosteriors.D <- t(ListSNPs$Marginals[[2]])
	ListSNPs.MarginalPosteriors.D <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors.D))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors.D[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors.D[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors.D) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.D.Melted <- melt(ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.D.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.D.Melted$Marker, levels=ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.I <- t(ListSNPs$Marginals[[3]])
	ListSNPs.MarginalPosteriors.I <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors.I))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors.I[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors.I[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors.I) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.I.Melted <- melt(ListSNPs.MarginalPosteriors.I[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.I.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.I.Melted$Marker, levels=ListSNPs.MarginalPosteriors.I[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.U <- t(ListSNPs$Marginals[[1]])
	ListSNPs.MarginalPosteriors.U <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors.U))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors.U[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors.U[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors.U) <- c("Marker", DataSources)	

	ListSNPs.MarginalPosteriors.U.Melted <- melt(ListSNPs.MarginalPosteriors.U[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.U.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.U.Melted$Marker, levels=ListSNPs.MarginalPosteriors.U[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	Z <- c()
	for (i in DataSources) {
		eval(parse(text=paste("Z <- cbind(Z, ListSNPs$SNPs$", i, "_ZScore)", sep="")))
	}

	Z.Directions <- apply(Z, c(1,2), function(x) { y <- 0; if (x > 0) { y <- 1; } else if (x < 0) { y <- -1; } else if (x == 0) { y <- 0; } else { y <- NA; }; return(y);})

	ListSNPs.MarginalPosteriors.D.wDirection <- ListSNPs.MarginalPosteriors.D
	ListSNPs.MarginalPosteriors.D.wDirection[,2:(length(DataSources)+1)] <- ListSNPs.MarginalPosteriors.D.wDirection[,2:(length(DataSources)+1)] * Z.Directions

	ListSNPs.MarginalPosteriors.D.wDirection.Melted <- melt(ListSNPs.MarginalPosteriors.D.wDirection[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.D.wDirection.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.D.wDirection.Melted$Marker, levels=ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.I.wDirection <- ListSNPs.MarginalPosteriors.I
	ListSNPs.MarginalPosteriors.I.wDirection[,2:(length(DataSources)+1)] <- ListSNPs.MarginalPosteriors.I.wDirection[,2:(length(DataSources)+1)] * Z.Directions

	ListSNPs.MarginalPosteriors.I.wDirection.Melted <- melt(ListSNPs.MarginalPosteriors.I.wDirection[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.I.wDirection.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.I.wDirection.Melted$Marker, levels=ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	ListSNPs.MarginalPosteriors.U.wDirection <- ListSNPs.MarginalPosteriors.U
	ListSNPs.MarginalPosteriors.U.wDirection[,2:(length(DataSources)+1)] <- ListSNPs.MarginalPosteriors.U.wDirection[,2:(length(DataSources)+1)] * Z.Directions

	ListSNPs.MarginalPosteriors.U.wDirection.Melted <- melt(ListSNPs.MarginalPosteriors.U.wDirection[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,])
	ListSNPs.MarginalPosteriors.U.wDirection.Melted$Marker <- factor(ListSNPs.MarginalPosteriors.U.wDirection.Melted$Marker, levels=ListSNPs.MarginalPosteriors.D[hclust(dist(ListSNPs.MarginalPosteriors.D[,2:(length(DataSources)+1)]))$order,1])

	if (PrintMarginals == 0) {
		if (Marginal == "U") {
			return(ggplot(ListSNPs.MarginalPosteriors.U.wDirection.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
		} else if (Marginal == "D") {
			return(ggplot(ListSNPs.MarginalPosteriors.D.wDirection.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
		} else if (Marginal == "I") {
			return(ggplot(ListSNPs.MarginalPosteriors.I.wDirection.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
		} else {
			return(PH <- 1)
		}
	} else if (PrintMarginals == 1) {
#		return(ListSNPs.MarginalPosteriors)
		return(PH <- 1)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}

}

#20170622 NOTE -- added/made this function below
PlotZScores.OrderedByClusteredMarginals.ForManuscript <- function (DataSources, ListSNPs, Marginal, PrintZScores=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="darkred", midcolor="white", highcolor="steelblue") {
	Marginal.Position <- -9
	
	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}
	
	ListSNPs.MarginalPosteriors <- t(ListSNPs$Marginals[[Marginal.Position]])
	ListSNPs.MarginalPosteriors <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors) <- c("Marker", DataSources)	

	ZScores <- data.frame(ListSNPs$SNPs$Marker)
	for (i in DataSources) {
		eval(parse(text=paste("ZScores <- cbind(ZScores, abs(ListSNPs$SNPs$", i, "_ZScore))", sep="")))
	}
	colnames(ZScores) <- c("Marker", DataSources)	

	ZScores.Melted <- melt(ZScores[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,])
	ZScores.Melted$Marker <- factor(ZScores.Melted$Marker, levels=ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,1])

	ZScores <- ZScores[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,]
	
	if (PrintZScores == 0) {
		return(ggplot(ZScores.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue") + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
	} else if (PrintZScores == 1) {
		return(ZScores)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}
}

#20170622 NOTE -- added/made this function below
PlotZScores.OrderedByClusteredMarginals.withDirection.ForManuscript <- function (DataSources, ListSNPs, Marginal, PrintZScores=1, Main=NULL, xLab=NULL, yLab=NULL, lowcolor="darkred", midcolor="white", highcolor="steelblue") {
	Marginal.Position <- -9
	
	if (Marginal == "U") {
		Marginal.Position <- 1
	} else if (Marginal == "D") {
		Marginal.Position <- 2
	} else if (Marginal == "I") {
		Marginal.Position <- 3
	} else {
		#20170219 CHECK_0 -- Prob: Put an exit of some sort here if Marginal not equal to any of these three options
		PH <- 1
	}
	
	ListSNPs.MarginalPosteriors <- t(ListSNPs$Marginals[[Marginal.Position]])
	ListSNPs.MarginalPosteriors <- data.frame(cbind(as.character(ListSNPs$SNPs$Marker), ListSNPs.MarginalPosteriors))
	for (i in 1:length(DataSources)) {
		ListSNPs.MarginalPosteriors[,i+1] <- as.numeric(as.character(ListSNPs.MarginalPosteriors[,i+1]))
	}
	colnames(ListSNPs.MarginalPosteriors) <- c("Marker", DataSources)	

	ZScores <- data.frame(ListSNPs$SNPs$Marker)
	for (i in DataSources) {
		eval(parse(text=paste("ZScores <- cbind(ZScores, ListSNPs$SNPs$", i, "_ZScore)", sep="")))
	}
	colnames(ZScores) <- c("Marker", DataSources)	

	ZScores.Melted <- melt(ZScores[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,])
	ZScores.Melted$Marker <- factor(ZScores.Melted$Marker, levels=ListSNPs.MarginalPosteriors[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,1])

	ZScores <- ZScores[hclust(dist(ListSNPs.MarginalPosteriors[,2:(length(DataSources)+1)]))$order,]
	
	if (PrintZScores == 0) {
		return(ggplot(ZScores.Melted, aes(variable, Marker, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title=Main))
	} else if (PrintZScores == 1) {
		return(ZScores)
	} else {
		#Print some type of error sign here?
		return(PH <- 1)
	}
}



