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

#Annotating merged data with, if provided, GWAS SNPs
#~~~~~~

AnnotateDataWithGWASSNPs_Old <- function (MergedDataSource1, GWASsnps1, BPWindow=500000) {
        GWASannot1 <- 0
#        print(MergedDataSource1)
	for (snpIndex in 1:nrow(GWASsnps1)) {
		if (GWASsnps1[snpIndex,]$Chr == as.numeric(as.character(MergedDataSource1["Chr"]))) {
                        if (GWASsnps1[snpIndex,]$BP == as.numeric(as.character(MergedDataSource1["BP"]))) {
                                GWASannot1 <- 1
                        }
                        else if ((GWASsnps1[snpIndex,]$BP >= as.numeric(as.character(MergedDataSource1["BP"])) - BPWindow) && (GWASsnps1[snpIndex,]$BP <= as.numeric(as.character(MergedDataSource1["BP"])) + BPWindow) && (GWASannot1 != 1)) {
                                GWASannot1 <- 2
                        }
                        else {
                                PH <- NULL
                        }
                }
        }
        return(GWASannot1)
}

AnnotateDataWithGWASSNPs <- function (MergedDataSource1, GWASsnps1, BPWindow=500000) {
        GWASannot1 <- rep(0, nrow(MergedDataSource1))
#        print(MergedDataSource1)
        for (snpIndex in 1:nrow(GWASsnps1)) {
        	GWASannot2 <- rep(0, nrow(MergedDataSource1))
        	GWASannot3 <- rep(0, nrow(MergedDataSource1))
		GWASannot2[MergedDataSource1$Chr == GWASsnps1[snpIndex,]$Chr & MergedDataSource1$BP - BPWindow <= GWASsnps1[snpIndex,]$BP & MergedDataSource1$BP + BPWindow >= GWASsnps1[snpIndex,]$BP] <- 2
		GWASannot2[MergedDataSource1$Chr == GWASsnps1[snpIndex,]$Chr & MergedDataSource1$BP == GWASsnps1[snpIndex,]$BP] <- 1
		
		GWASannot3[GWASannot1==2 | GWASannot2 == 2] <- 2
		GWASannot3[GWASannot1==1 | GWASannot2 == 1] <- 1

		GWASannot1 <- GWASannot3
	}
        return(GWASannot1)
}

GetZScoreAndDirection <- function(DataSources1) {
        ZScore <- qnorm(log(as.numeric(as.character(DataSources1["pValue"]))/2), lower.tail=FALSE, log.p=TRUE);
        if (DataSources1["Direction"] == "-") {
                ZScore <- ZScore * -1;
        }
        return(ZScore);
}

CheckForInfiniteZScores <- function(DataSources1_ZScores) {
        returnValue <- FALSE
        CheckValues <- is.infinite(DataSources1_ZScores)
        if (TRUE %in% CheckValues) {
                returnValue <- TRUE
        }
        return(returnValue)
}

ReplaceInfiniteZScoresWithMax <- function(DataSources1_ZScores) {
        CheckValues <- is.infinite(DataSources1_ZScores)
        if (TRUE %in% CheckValues) {
                CheckValues_positive <- is.infinite(DataSources1_ZScores) & DataSources1_ZScores > 0
                CheckValues_negative <- is.infinite(DataSources1_ZScores) & DataSources1_ZScores < 0
                maxZScore <- max(DataSources1_ZScores[!CheckValues])
                minZScore <- min(DataSources1_ZScores[!CheckValues])
                DataSources1_ZScores[CheckValues_positive] <- maxZScore
                DataSources1_ZScores[CheckValues_negative] <- minZScore
        }
        return(DataSources1_ZScores)
}

MergeDataSources <- function (DataSources, LogFile) {

	#Preparing and merging data
        #~~~~~~

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Beginning DataSources merging.", sep=""))

        MergedDataSources <- data.frame()
        for (CurrentDataSource in DataSources) {

                #20160930 CHECK_0 -- Prob: Go over and do this/figure out what want to do?
		###### do thissssss
                #
                #Keep MAF and average across datasets? Check A1 for consistency or not expecting that since differing directions?
                #
                ###### do thissssss

                if (nrow(MergedDataSources)==0) {
                        MergedDataSources <- eval(parse(text=CurrentDataSource))
                        MergedDataSources$ZScore <- apply(MergedDataSources[,c("pValue", "Direction")], 1, GetZScoreAndDirection)
#                       MergedDataSources$pValue <- NULL
#                       MergedDataSources$Direction <- NULL
                        MergedDataSources_namesCurrent <- names(MergedDataSources)
                        MergedDataSources_namesNew <- c()
                        for (columnHeader1 in MergedDataSources_namesCurrent) {
                                if (columnHeader1 %in% c("Chr", "BP", "A1", "MAF")) {
                                        MergedDataSources_namesNew <- c(MergedDataSources_namesNew, columnHeader1)
                                }
                                else {
                                        MergedDataSources_namesNew <- c(MergedDataSources_namesNew, paste(CurrentDataSource, "_", columnHeader1, sep=""))
                                }
                        }
                        names(MergedDataSources) <- MergedDataSources_namesNew
                        MergedDataSources$ChrBP <- paste(MergedDataSources$Chr, MergedDataSources$BP, sep="_")
                }
                else {
                        CurrentDataSource_temp <- eval(parse(text=paste(CurrentDataSource, "[,c(\"Chr\", \"BP\", \"Direction\", \"pValue\", \"N\")]", sep="")))
                        CurrentDataSource_temp$ZScore <- apply(CurrentDataSource_temp[,c("pValue", "Direction")], 1, GetZScoreAndDirection)
#                       CurrentDataSource_temp$Direction <- NULL
#                       CurrentDataSource_temp$pValue <- NULL
                        CurrentDataSource_temp_namesCurrent <- names(CurrentDataSource_temp)
                        CurrentDataSource_temp_namesNew <- c()
                        for (columnHeader1 in CurrentDataSource_temp_namesCurrent) {
                                CurrentDataSource_temp_namesNew <- c(CurrentDataSource_temp_namesNew, paste(CurrentDataSource, "_", columnHeader1, sep=""))
                        }
                        names(CurrentDataSource_temp) <- CurrentDataSource_temp_namesNew
                        CurrentDataSource_temp$ChrBP <- paste(eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_Chr", sep=""))), eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_BP", sep=""))), sep="_")
                        eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_Chr <- NULL", sep="")))
                        eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_BP <- NULL", sep="")))

                	#20160930 CHECK_0 -- Prob: Do this/figure this out? If even wanting to do something like this?
                        ###### do thissssss
                        #
                        # Test performances of merge function here in varying circumstances
                        #
                        ###### do thissssss

			###Candidate for unit tests
                        MergedDataSources <- merge(MergedDataSources, CurrentDataSource_temp, by="ChrBP")

                        rm(CurrentDataSource_temp)
                }
        }

	#Checking all SNPs are now oriented to minor allele, and making associated changes where needed (eg flipping MAF, ZScore directions)
	#20161108 CHECK_0 -- Prob: Ask users to give both A1 and A2 so that can flip properly if needed here
	MAF_CheckList <- MergedDataSources$MAF > .5
	MergedDataSources[MAF_CheckList,]$MAF <- 1 - MergedDataSources[MAF_CheckList,]$MAF
	for (CurrentDataSource in DataSources) {
		eval(parse(text=paste("MergedDataSources[MAF_CheckList,]$", CurrentDataSource, "_ZScore <- -1 * MergedDataSources[MAF_CheckList,]$", CurrentDataSource, "_ZScore", sep="")))
	}

	return(list(MergedDataSources=MergedDataSources, LogFile=LogFile))

}

AnnotateMergedDataWithGWASSNPs <- function(MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, LogFile) {
        
	#Annotating merged data with, if provided, GWAS SNPs
        #~~~~~~

        if (is.null(GWASsnps)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- No GWASsnps list provided, skipping annotating MergedDataSources.", sep=""))
                MergedDataSources$GWASannot <- 0
        }
        else {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Annotating MergedDataSources with provided GWASsnps list.", sep=""))
#                MergedDataSources$GWASannot <- apply(MergedDataSources, 1, AnnotateDataWithGWASSNPs_Old, GWASsnps1=GWASsnps, BPWindow=GWASsnps_AnnotateWindow)
                MergedDataSources$GWASannot <- AnnotateDataWithGWASSNPs(MergedDataSources[,c("Chr", "BP")], GWASsnps, GWASsnps_AnnotateWindow)
        }

	return(list(MergedDataSources=MergedDataSources, LogFile=LogFile))

}

ProcessMergedAndAnnotatedDataSources <- function (DataSources, MergedDataSources, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, LogFile) {

        #Calculating RSS0(eg ZScoreCorMatrix) and subsetting down to marginally significant SNPs
        #~~~~~~
        
	MarginalSNPs <- list()

        #Getting ZScore matrix
        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Determining Z-score correlation matrix.", sep=""))

        ZScoresFull_CommandText <- ""
        for (DataSource in DataSources) {
                if (length(strsplit(ZScoresFull_CommandText, "")[[1]]) == 0) {
                        ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, "cbind(MergedDataSources$", DataSource, "_ZScore", sep="")
                }
                else {
                        ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, ",MergedDataSources$", DataSource, "_ZScore", sep="")
                }
        }
        ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, ")", sep="")
        ZScoresFull <- eval(parse(text=ZScoresFull_CommandText))

        ZScoresFullNames_CommandText <- c()
        for (DataSource in DataSources) {
                ZScoresFullNames_CommandText <- c(ZScoresFullNames_CommandText, paste(DataSource, "_ZScore", sep=""))
        }
        colnames(ZScoresFull) <- ZScoresFullNames_CommandText
#       print(ZScoresFull)

        #Checking ZScore matrix for infinites and replacing with appropriate max/min values
        ZScoresFull_InfiniteCheck <- apply(ZScoresFull, 2, CheckForInfiniteZScores)
        if (TRUE %in% ZScoresFull_InfiniteCheck) {
                warning(paste(format(Sys.time()), " -- One of your datafiles has p-values less than the threshold which R can properly convert them to log-scale, meaning their log-values return as 'Infinite'. bmass automatically replaces these 'Infinite' values with the max, non-infinite Z-scores available, but it is recommended you self-check this as well.", sep=""))
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- One of your datafiles has p-values less than the threshold which R can properly convert them to log-scale, meaning their log-values return as 'Infinite'. bmass automatically replaces these 'Infinite' values with the max, non-infinite Z-scores available, but it is recommended you self-check this as well.", sep=""))
                ZScoresFull <- apply(ZScoresFull, 2, ReplaceInfiniteZScoresWithMax)
                for (ZScoresFull_colName in colnames(ZScoresFull)) {
                        eval(parse(text=paste("MergedDataSources$", ZScoresFull_colName, " <- ZScoresFull[,\"", ZScoresFull_colName, "\"]", sep="")))
                }
        }

        #Getting ZScore null subset matrix, calculating correlation matrix, and associated naive multivariate statistics
        ZScoresNullSetSelection_CommandText <- ""
        for (DataSource in DataSources) {
                if (length(strsplit(ZScoresNullSetSelection_CommandText, "")[[1]]) == 0) {
                        ZScoresNullSetSelection_CommandText <- paste(ZScoresNullSetSelection_CommandText, "(abs(ZScoresFull[,\"", DataSource, "_ZScore\"])<2)", sep="")
                }
                else {
                        ZScoresNullSetSelection_CommandText <- paste(ZScoresNullSetSelection_CommandText, " & (abs(ZScoresFull[,\"", DataSource, "_ZScore\"])<2)", sep="")
                }
        }
        ZScoresNullSetSelection <- eval(parse(text=ZScoresNullSetSelection_CommandText))
        ZScoresNullset <- ZScoresFull[ZScoresNullSetSelection,]

        ZScoresCorMatrix <- cor(ZScoresNullset)

        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Determining initial threshold statistics.", sep=""))

        ZScoresCorMatrix_Inverse <- chol2inv(chol(ZScoresCorMatrix))
        ZScoresFull_mvstat <- rowSums(ZScoresFull * (ZScoresFull %*% ZScoresCorMatrix_Inverse))
        MergedDataSources$mvstat <- ZScoresFull_mvstat
        ZScoresFull_mvstat_log10pVal <- -log10(exp(1))*pchisq(ZScoresFull_mvstat, df=ncol(ZScoresFull), log.p=TRUE, lower.tail=FALSE)
        MergedDataSources$mvstat_log10pVal <- ZScoresFull_mvstat_log10pVal

        ZScoresFull_unistat <- apply(ZScoresFull^2, 1, max)
        MergedDataSources$unistat <- ZScoresFull_unistat
        ZScoresFull_unistat_log10pVal <- -log10(exp(1))*pchisq(ZScoresFull_unistat, df=1, log.p=TRUE, lower.tail=FALSE)
        MergedDataSources$unistat_log10pVal <- ZScoresFull_unistat_log10pVal

#       print(MergedDataSources)

        #Creating subset of marginally significant SNPs using SNPMarginalUnivariateThreshold and SNPMarginalMultivariateThreshold
        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Subsetting down to marginally significant SNPs based on univariate and multivariate thresholds: ", as.character(SNPMarginalUnivariateThreshold) ," & ", as.character(SNPMarginalMultivariateThreshold) ,".", sep=""))
	MarginalSNPs$SNPs <- MergedDataSources[MergedDataSources$mvstat_log10pVal > -log10(SNPMarginalUnivariateThreshold) | MergedDataSources$unistat_log10pVal > -log10(SNPMarginalMultivariateThreshold),]

#       print(MarginalSNPs)
	
	return(list(MarginalSNPs=MarginalSNPs, ZScoresCorMatrix=ZScoresCorMatrix, LogFile=LogFile))

}



