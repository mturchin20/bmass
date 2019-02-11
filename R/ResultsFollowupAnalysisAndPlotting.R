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

#' @title Get Model Prior Matrix
#'
#' @description Creates a matrix containing the model descriptions and
#' their associated priors.
#' 
#' @param DataSources A string indicating the variable names of the
#' input datafiles and phenotypes.
#' 
#' @param Models A matrix describing the models being explored
#' (default output from running bmass).
#' 
#' @param ModelPriors A vector containing the priors on each model
#' across each tranche of sigma alpha (default output from running
#' bmass; length is number of models times number of of sigma alphas).
#' 
#' @param LogFile A matrix of string outputs for function logging
#' purposes (default output from running bmass).
#' 
#' @param SigmaAlphas A vector containing the different values
#' traversed for this 'effect size controlling' hyperparameter (see
#' "Prior on Sigma_Alpha" in Stephens 2013 PLoS ONE.
#'
#' @return A matrix containing the original description of each model
#' sort by prior, each model's trained prior, the cummulative prior
#' distribution, and the model's original order position.
#"
#' @examples
#' \dontrun{
#' GetModelPriorMatrix(c("HDL", "LDL", "TG", "TC"), bmassOutput$Models, bmassOutput$ModelPriors, bmassOutput$LogFile)
#' bmassOutput[c("ModelPriorMatrix", "LogFile")] <- GetModelPriorMatrix(c("HDL", "LDL", "TG", "TC"), bmassOutput$Models, bmassOutput$ModelPriors, bmassOutput$LogFile)
#' }
#'
#' @export
#' 
GetModelPriorMatrix <- function (DataSources, Models, ModelPriors, LogFile, SigmaAlphas = c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)) {

	ModelPriors.Collapsed <- CollapseSigmaAlphasTogether(ModelPriors, length(SigmaAlphas))

	ModelPriorMatrix <- cbind(Models,ModelPriors.Collapsed)[order(ModelPriors.Collapsed, decreasing=TRUE),]
	ModelPriorMatrix <- data.frame(cbind(ModelPriorMatrix, cumsum(ModelPriorMatrix[,(length(DataSources)+1)])), order(ModelPriors.Collapsed, decreasing=TRUE))
	colnames(ModelPriorMatrix) <- c(DataSources, "Prior", "Cumm_Prior", "OrigOrder")

	return(list(ModelPriorMatrix=ModelPriorMatrix, LogFile=LogFile))

}

#' @title Get Top Multivariate Models
#'
#' @description Get a summary of the top models per SNP across all
#' multivariate \{U,D,I\} combinations based on posterior
#' probabilities.
#' 
#' @param DataSources A string indicating the variable names of the
#' input datafiles and phenotypes.
#' 
#' @param ListSNPs A list produced from running bmass containing the
#' SNPs of interest to get marginal posteriors for.
#' 
#' @param ModelPriorMatrix A matrix detailing the models being
#' explored and their associated priors (obtained by running
#' 'GetModelPriorMatrix()')
#' 
#' @param LogFile A matrix of string outputs for function logging
#' purposes (default output from running bmass).
#'
#' @return A matrix containing each model that was a SNP's top model
#' at least once, along with related information; this matrix is
#' appended to the input ListSNPs as a new object, "TopModels" (the
#' full returned object is a list containing the input ListSNPs and
#' the input LogFile).
#'
#' @examples
#' \dontrun{
#' GetTopModelsPerSNPViaPosteriors(c("HDL", "LDL", "TG", "TC"), bmassOutput$NewSNPs, bmassOutput$ModelPriorMatrix, bmassOutput$LogFile)
#' bmassOutput[c("NewSNPs", "LogFile")] <- GetTopModelsPerSNPViaPosteriors(c("HDL", "LDL", "TG", "TC"), bmassOutput$NewSNPs, bmassOutput$ModelPriorMatrix, bmassOutput$LogFile)
#' }
#'
#' @export
#' 
GetTopModelsPerSNPViaPosteriors <- function (DataSources, ListSNPs, ModelPriorMatrix, LogFile) {

	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Running GetTopModelsPerSNPViaPosteriors().", sep=""))
	
	ListSNPs.Posteriors <- cbind(ListSNPs$Posteriors[,(length(DataSources)+1):ncol(ListSNPs$Posteriors)])[ModelPriorMatrix[,ncol(ModelPriorMatrix)],]
	
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

#' @title Get Marginal \{U,D,I\} Posteriors
#'
#' @description Get marginal posteriors for how much every individual
#' phenotype belongs to categories \{U,D,I\} across each SNP
#' 
#' @param DataSources A string indicating the variable names of the
#' input datafiles and phenotypes.
#' 
#' @param ListSNPs A list produced from running bmass containing the
#' SNPs of interest to get marginal posteriors for.
#' 
#' @param Models A matrix describing the models being explored
#' (default output from running bmass).
#' 
#' @param LogFile A matrix of string outputs for function logging
#' purposes (default output from running bmass).
#'
#' @return A list containing three matrices of SNPs x Phenotypes
#' marginal posteriors for each category \{U,D,I\}; this list is
#' appended to the input ListSNPs as a new object, "Posteriors" (the
#' full returned object is a list containing the input ListSNPs and
#' the input LogFile).
#'
#' @examples
#' \dontrun{
#' GetMarginalPosteriors(c("HDL", "LDL", "TG", "TC"), bmassOutput$NewSNPs, bmassOutput$Models, bmassOutput$LogFile)
#' bmassOutput[c("NewSNPs", "LogFile")] <- GetMarginalPosteriors(c("HDL", "LDL", "TG", "TC"), bmassOutput$NewSNPs, bmassOutput$Models, bmassOutput$LogFile)
#' }
#'
#' @export
#' 
GetMarginalPosteriors <- function (DataSources, ListSNPs, Models, LogFile) {
	
	LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Running GetMarginalPosteriors().", sep=""))
        
	ListSNPs.MarginalPosteriors <- marginal.postprobs(ListSNPs$Posteriors[,(length(DataSources)+1):ncol(ListSNPs$Posteriors)], Models, 1)	
	rownames(ListSNPs.MarginalPosteriors$pU) <- DataSources 

	ListSNPs$Marginals <- ListSNPs.MarginalPosteriors

	return(list(ListSNPs=ListSNPs, LogFile=LogFile))

}
