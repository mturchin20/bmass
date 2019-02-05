AnnotateDataWithGWASSNPs <- function (MergedDataSource1, GWASsnps1, BPWindow=500000) {
        GWASannot1 <- rep(0, nrow(MergedDataSource1))
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

#' GetZScoreAndDirection
#'
#' GetZScoreAndDirection
#' 
#' @keywords internal
#'
#' @importFrom stats qnorm
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

                if (nrow(MergedDataSources)==0) {
                        MergedDataSources <- eval(parse(text=CurrentDataSource))
                        MergedDataSources$ZScore <- apply(MergedDataSources[,c("pValue", "Direction")], 1, GetZScoreAndDirection)
#                       MergedDataSources$pValue <- NULL
#                       MergedDataSources$Direction <- NULL
                        MergedDataSources_namesCurrent <- names(MergedDataSources)
                        MergedDataSources_namesNew <- c()
                        for (columnHeader1 in MergedDataSources_namesCurrent) {
                                if (columnHeader1 %in% c("Chr", "BP", "Marker", "A1", "MAF")) {
                                        MergedDataSources_namesNew <- c(MergedDataSources_namesNew, columnHeader1)
				} else {
                                        MergedDataSources_namesNew <- c(MergedDataSources_namesNew, paste(CurrentDataSource, "_", columnHeader1, sep=""))
                                }
                        }
                        names(MergedDataSources) <- MergedDataSources_namesNew
                        MergedDataSources$ChrBP <- paste(MergedDataSources$Chr, MergedDataSources$BP, sep="_")
		} else {
                        CurrentDataSource_temp <- eval(parse(text=paste(CurrentDataSource, "[,c(\"Chr\", \"BP\", \"A1\", \"Direction\", \"pValue\", \"N\")]", sep="")))
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
			
#			if (FALSE %in% eval(parse(text=paste("ifelse(MergedDataSources$A1==MergedDataSources$", CurrentDataSource, "_A1), TRUE, FALSE)))) {
#
#			}

			eval(parse(text=paste("CurrentDataSource_temp$", CurrentDataSource, "_A1 <- NULL", sep="")))
                        MergedDataSources <- merge(MergedDataSources, CurrentDataSource_temp, by="ChrBP")


                        rm(CurrentDataSource_temp)
                }
        }

	#Checking that SNP remains MAF > .5 after processing
	MAF_CheckList <- MergedDataSources$MAF > .5
	if (TRUE %in% MAF_CheckList) {
		stop(Sys.time(), " -- the created MergedDataSources file has SNPs with final MAFs > .5. Please fix and rerun bmass: ", MergedDataSources[MAF_CheckList,]$Marker) 
	}

	return(list(MergedDataSources=MergedDataSources, LogFile=LogFile))

}

AnnotateMergedDataWithGWASSNPs <- function(MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, LogFile) {
        
	#Annotating merged data with, if provided, GWAS SNPs
        #~~~~~~

        if (is.null(GWASsnps)) {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- No GWASsnps list provided, skipping annotating MergedDataSources.", sep=""))
                MergedDataSources$GWASannot <- 0
        } else {
                LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Annotating MergedDataSources with provided GWASsnps list.", sep=""))
                MergedDataSources$GWASannot <- AnnotateDataWithGWASSNPs(MergedDataSources[,c("Chr", "BP")], GWASsnps, GWASsnps_AnnotateWindow)
        }

	return(list(MergedDataSources=MergedDataSources, LogFile=LogFile))

}

#' ProcessMergedAndAnnotatedDataSources
#'
#' ProcessMergedAndAnnotatedDataSources
#' 
#' @keywords internal
#'
#' @importFrom stats cor pchisq
ProcessMergedAndAnnotatedDataSources <- function (DataSources, MergedDataSources, ZScoresCorMatrix, SNPMarginalUnivariateThreshold, SNPMarginalMultivariateThreshold, LogFile) {

        #Calculating RSS0(eg ZScoreCorMatrix) and subsetting down to marginally significant SNPs
        #~~~~~~
        
	MarginalSNPs <- list()

        #Getting ZScore matrix
        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Determining Z-score correlation matrix.", sep=""))

        ZScoresFull_CommandText <- ""
        for (DataSource in DataSources) {
                if (length(strsplit(ZScoresFull_CommandText, "")[[1]]) == 0) {
                        ZScoresFull_CommandText <- paste(ZScoresFull_CommandText, "cbind(MergedDataSources$", DataSource, "_ZScore", sep="")
                } else {
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
                } else {
                        ZScoresNullSetSelection_CommandText <- paste(ZScoresNullSetSelection_CommandText, " & (abs(ZScoresFull[,\"", DataSource, "_ZScore\"])<2)", sep="")
                }
        }
        ZScoresNullSetSelection <- eval(parse(text=ZScoresNullSetSelection_CommandText))
        ZScoresNullset <- ZScoresFull[ZScoresNullSetSelection,]

	#If no input `ZScoresCorMatrix`, then should be `NULL`
        if (is.null(ZScoresCorMatrix)) {
		ZScoresCorMatrix <- cor(ZScoresNullset)
	}

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

        #Creating subset of marginally significant SNPs using SNPMarginalUnivariateThreshold and SNPMarginalMultivariateThreshold
        LogFile <- rbind(LogFile, paste(format(Sys.time()), " -- Subsetting down to marginally significant SNPs based on univariate and multivariate thresholds: ", as.character(SNPMarginalUnivariateThreshold) ," & ", as.character(SNPMarginalMultivariateThreshold) ,".", sep=""))
	MarginalSNPs$SNPs <- MergedDataSources[MergedDataSources$mvstat_log10pVal > -log10(SNPMarginalUnivariateThreshold) | MergedDataSources$unistat_log10pVal > -log10(SNPMarginalMultivariateThreshold),]

	return(list(MergedDataSources=MergedDataSources, MarginalSNPs=MarginalSNPs, ZScoresCorMatrix=ZScoresCorMatrix, LogFile=LogFile))

}



