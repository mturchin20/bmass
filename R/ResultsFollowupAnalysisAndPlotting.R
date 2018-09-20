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

GetModelPriorMatrix <- function (DataSources, Models, ModelPriors, SigmaAlphas, LogFile) {

	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- stuff1.", sep=""))

	ModelPriors.Collapsed <- CollapseSigmaAlphasTogether(ModelPriors, length(SigmaAlphas))

	ModelPriorMatrix <- cbind(Models,ModelPriors.Collapsed)[order(ModelPriors.Collapsed, decreasing=TRUE),]
	ModelPriorMatrix <- data.frame(cbind(ModelPriorMatrix, cumsum(ModelPriorMatrix[,(length(DataSources)+1)])), order(ModelPriors.Collapsed, decreasing=TRUE))
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

	#The below code leads to an error downstream when there is only 1 entry, hence the if clause and separation (code is also unnecessary if there is only 1 entry)
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

