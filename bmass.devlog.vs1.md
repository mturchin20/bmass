##General Comments

###Markdown related

Including spaces/linebreaks: http://stackoverflow.com/questions/24575680/new-lines-inside-paragraph-in-readme-md

##20160812

Copying in original code from Matthew's multivariate git repo (https://github.com/stephens999/multivariate) with names that signify their original source. Going to use them as a starting point for editing/cleaning the main function files

```
cp -p /Users/mturchin20/Documents/Work/LabMisc/StephensLab/Multivariate/multivariate/test.funcs.R /Users/mturchin20/Documents/Work/LabMisc/StephensLab/bmass/R/Stephens2013PLoSONE.test.funcs.R
cp -p /Users/mturchin20/Documents/Work/LabMisc/StephensLab/Multivariate/multivariate/globallipids/GLCfuncs.R /Users/mturchin20/Documents/Work/LabMisc/StephensLab/bmass/R/Stephens2013PLoSONE.globallipids.GLCfuncs.R
```


##20161006

###Including tags onto current and past commits

```
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git tag -a v0.1.1 -m "First draft of code reorganization, bmass() runs"
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git push origin master
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git log --pretty=oneline | head -n 10
4134d2dc370ad884f7ac597ccd41aa35fde9d41a 20161006 -- Updating NEWS.md
b00064673d257dec6c208187f100256e205805b1 20161006 -- Pushing v0.1.1: First draft of reorganized code & the bmass() function runs
01b70ea42ca1c3dd4bbcb811007d02f4d1986605 20161006 -- Updating some main .R files
6fb16ec275a6de5182856512cc5d751d75a08393 20161006 -- Updating some main .R files
c8a904535486c9e495125be65fc0083b208389a8 20161006 -- Updating main .R files mid reorganization efforts
8fe260e661fff1adf9199df312ccdd48b8e50ad9 20161005 -- Updating main .R files for reorganization efforts
180000f75a0d0c416d68e8662a4d91f2b9f2bf9a 20160930 -- Pushing updates for unit tests as well as continued reorganization of code chunks
d7dcc330858aea380e837990ade11cace6820e85 20160930 -- Pushing updates for unit tests, first complete revision for test-CheckDataSources.R
2a11a4ee39df21ef97477268d05d8492a72db967 20160928 -- Pushing first few complete unit tests; they successfully run via devtools::test()
38edfbd787af0f48c798ac46a99fc3ff5451fe17 20160928 -- checking in updates as begin work on unit tests, test datasets, and associated Roxygen-based documentation
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git log --pretty=oneline | vi -
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git tag -a v0.1.0 -m "First succesful run of bmass()" fcc9e63df5f5f4dda182bea7a5888ce66556d7bf
```

Relevant URLs to the above:
<br />
https://git-scm.com/book/en/v2/Git-Basics-Tagging  
http://stackoverflow.com/questions/18216991/create-a-tag-in-github-repository  
http://stackoverflow.com/questions/5480258/how-to-delete-a-remote-tag  

<br />
###Practicing/working on R CMD BUILD/CHECK functions

```
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab]$R CMD build bmass
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab]$R CMD CHECK bmass_0.1.1.9000.tar.gz 
* using log directory '/Users/mturchin20/Documents/Work/LabMisc/StephensLab/bmass.Rcheck'
* using R version 3.3.0 (2016-05-03)
* using platform: x86_64-apple-darwin15.4.0 (64-bit)
* using session charset: ASCII
* checking for file 'bmass/DESCRIPTION' ... OK
* this is package 'bmass' version '0.1.1.9000'
* package encoding: UTF-8
* checking package namespace information ... OK
* checking package dependencies ... OK
* checking if this is a source package ... OK
* checking if there is a namespace ... OK
* checking for executable files ... OK
* checking for hidden files and directories ... OK
* checking for portable file names ... OK
* checking for sufficient/correct file permissions ... OK
* checking whether package 'bmass' can be installed ... OK
* checking installed package size ... OK
* checking package directory ... OK
* checking DESCRIPTION meta-information ... OK
* checking top-level files ... OK
* checking for left-over files ... OK
* checking index information ... OK
* checking package subdirectories ... OK
* checking R files for non-ASCII characters ... OK
* checking R files for syntax errors ... OK
* checking whether the package can be loaded ... OK
* checking whether the package can be loaded with stated dependencies ... OK
* checking whether the package can be unloaded cleanly ... OK
* checking whether the namespace can be loaded with stated dependencies ... OK
* checking whether the namespace can be unloaded cleanly ... OK
* checking loading without being on the library search path ... OK
* checking dependencies in R code ... OK
* checking S3 generic/method consistency ... OK
* checking replacement functions ... OK
* checking foreign function calls ... OK
* checking R code for possible problems ... NOTE
computeprior: warning in tabulate(z + 1, nbin = 3): partial argument
  match of 'nbin' to 'nbins'
DetermineAndApplyPriors: no visible global function definition for
  'runif'
ExploreBestModelsUsingPosteriors: no visible binding for global
  variable 'Prior_PreviousSNPsEB'
ExploreBestModelsUsingPosteriors: no visible binding for global
  variable 'SigmaAlphas'
ExploreBestModelsUsingPosteriors: no visible binding for global
  variable 'Prior_PreviousSNPsEB_check2'
ExploreBestModelsUsingPosteriors: no visible binding for global
  variable 'PreviousSNPs_logBFs_Stacked'
ExploreBestModelsUsingPosteriors: no visible binding for global
  variable 'DataSources'
GetZScoreAndDirection: no visible global function definition for
  'qnorm'
ProcessMergedAndAnnotatedDataSources: no visible global function
  definition for 'cor'
ProcessMergedAndAnnotatedDataSources: no visible global function
  definition for 'pchisq'
WriteTableLogFile: no visible global function definition for
  'write.table'
lbf.all: no visible binding for global variable 'allones'
Undefined global functions or variables:
  DataSources PreviousSNPs_logBFs_Stacked Prior_PreviousSNPsEB
  Prior_PreviousSNPsEB_check2 SigmaAlphas allones cor pchisq qnorm
  runif write.table
Consider adding
  importFrom("stats", "cor", "pchisq", "qnorm", "runif")
  importFrom("utils", "write.table")
to your NAMESPACE file.
* checking Rd files ... OK
* checking Rd metadata ... OK
* checking Rd cross-references ... OK
* checking for missing documentation entries ... WARNING
Undocumented code objects:
  'TestData1.DoesNotThrowWarningsVersion'
  'TestData1.ThrowsWarningsVersion'
  'TestData2.DoesNotThrowWarningsVersion'
  'TestData2.ThrowsWarningsVersion' 'bmass_TestData1' 'bmass_TestData2'
  'bmass_TestSigSNPs'
Undocumented data sets:
  'TestData1.DoesNotThrowWarningsVersion'
  'TestData1.ThrowsWarningsVersion'
  'TestData2.DoesNotThrowWarningsVersion'
  'TestData2.ThrowsWarningsVersion' 'bmass_TestData1' 'bmass_TestData2'
  'bmass_TestSigSNPs'
All user-level objects in a package should have documentation entries.
See chapter 'Writing R documentation files' in the 'Writing R
Extensions' manual.
* checking for code/documentation mismatches ... WARNING
Functions or methods with usage in documentation object 'CheckCharacterFormat' but not in code:
  CheckCharacterFormat

* checking Rd \usage sections ... WARNING
Undocumented arguments in documentation object 'CheckCharacterFormat'
  'DataSource1'
Documented arguments not in \usage in documentation object 'CheckCharacterFormat':
  'x' 'y' 'DataFileList' 'DataFileLocations' 'ExpectedColumnNames'

Undocumented arguments in documentation object 'OpenLogFile'
  'OutputFileBase'
Documented arguments not in \usage in documentation object 'OpenLogFile':
  'x' 'y'

Functions with \usage entries need to have the appropriate \alias
entries, and all their arguments documented.
The \usage entries must correspond to syntactically valid R code.
See chapter 'Writing R documentation files' in the 'Writing R
Extensions' manual.
* checking Rd contents ... OK
* checking for unstated dependencies in examples ... OK
* checking contents of 'data' directory ... OK
* checking data for non-ASCII characters ... OK
* checking data for ASCII and uncompressed saves ... OK
* checking examples ... WARNING
checking a package with encoding  'UTF-8'  in an ASCII locale

 ERROR
Running examples in 'bmass-Ex.R' failed
The error most likely occurred in:

> ### Name: CheckCharacterFormat
> ### Title: Checking and preparing input datafiles.
> ### Aliases: CheckCharacterFormat
> 
> ### ** Examples
> 
> func(1, 1)
Error: could not find function "func"
Execution halted
* checking for unstated dependencies in 'tests' ... OK
* checking tests ...
  Running 'testthat.R'
 OK
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
* checking PDF version of manual without hyperrefs or index ... ERROR
Re-running with no redirection of stdout/stderr.
Hmm ... looks like a package
Error in texi2dvi(file = file, pdf = TRUE, clean = clean, quiet = quiet,  : 
  pdflatex is not available
Error in texi2dvi(file = file, pdf = TRUE, clean = clean, quiet = quiet,  : 
  pdflatex is not available
Error in running tools::texi2pdf()
You may want to clean up by 'rm -rf /var/folders/hw/3fshyx49303c6gkw95j44pzh0000gn/T//Rtmpw2mUZG/Rd2pdf1ff4338c9296'
* DONE

Status: 2 ERRORs, 5 WARNINGs, 1 NOTE
See
  '/Users/mturchin20/Documents/Work/LabMisc/StephensLab/bmass.Rcheck/00check.log'
for details.
```

####Fixing ``no visible global function definition for 'runif'` 
From http://stackoverflow.com/questions/31132552/no-visible-global-function-definition-for-median, adding:

```
Imports:
    stats
```

to `DESCRIPTION` and `import(stats)` to NAMESPACE fixed this

This approach was also use to fix `no visible global function definition for 'write.table'`, where `utils` needed to be included in lieu of `stats`

####Other remaining misc for the time being

Most of the other comments in the `Note` are dealing with the not-completed function `ExploreBestModelsUsingPosteriors`, so a handful of variables are currently not initialized.  
Additionally, incomplete Roxygen comments are leading to one of the `Error`s, and the other `Error` is a product of pdflatex not being available? For the moment going to come back to finish these up later, but mostly look fine for now 


<br />
###Changed wrong commit author, following instructions from https://help.github.com/articles/changing-author-info/

TempScript.sh (copy/pasted from above with info below included):
```
OLD_EMAIL="mturchin20@uchicago.edu"
CORRECT_NAME="mturchin20"
CORRECT_EMAIL="mturchin20@gmail.com"
```

Command line:
```
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git clone https://github.com/mturchin20/bmass/
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git log | vi -
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$bash TempScript.sh
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git log | vi -
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git push --force --tags origin 'refs/heads/*'
##Then on other git repo location
#[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git pull origin master
#NOTE -- the above line was actually wrong, it fetched then merged which was a problem. As a result I had double the commit messages, duplicates of every one except for the one I changed. The below line helped me revert back to where I needed to start from again
#[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git reset --hard master@{"10 minutes ago"}
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git fetch --all
[  mturchin20@wireless-s0-no-150-8-172  ~/Documents/Work/LabMisc/StephensLab/bmass]$git reset --hard origin/master
```

Other websites used to help troubleshoot this:
<br />
http://stackoverflow.com/questions/1125968/how-to-force-git-pull-to-overwrite-local-files
http://stackoverflow.com/questions/1223354/undo-git-pull-how-to-bring-repos-to-old-state


<br />
<br />

##20161106

###Work related to moving to quicker version of `AnnotateDataWithGWASSNPs` 

Some initial exploratory runs of `AnnotateMergedDataWithGWASSNPs()` (which mainly calls/uses `AnnotateDataWithGWASSNPs`) to see how long it takes and what variables/factors influence it most. Appears, unsurprisingly, that there's a linear relationship between run time and number of merged dataset entries that need to be annotated.

```
> bmassOutput[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]
> system.time(bmassOutput[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
     user    system   elapsed
46465.448     2.272 46480.563
> dim(GWASsnps)
[1] 101   2
> GWASsnps_top25 <- head(GWASsnps, n=25)
> dim(GWASsnps_top25)
[1] 25  2
> system.time(bmassOutput[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps_top25, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
> bmassOutput2 <- bmassOutput
> system.time(bmassOutput2[c("MergedDataSources", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources, GWASsnps_top25, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")]
) 
     user    system   elapsed
11642.902     1.026 11646.142
> posterior.gl.prior <- posteriorprob(lbf.bigmat, normalize(rep(c(0,lbf$prior[-1]),nsigma)))
Error in t(lbf) : object 'lbf.bigmat' not found
```

```
> bmassOutput$MergedDataSources_Short <- bmassOutput$MergedDataSources[1:10000,]
.
.
.
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources_Short, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
   user  system elapsed
109.301   0.006 109.315
> bmassOutput$MergedDataSources_Short <- bmassOutput$MergedDataSources[1:100000,]                                                                                                                                              > system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources_Short, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
    user   system  elapsed
1072.733    0.076 1072.899
> bmassOutput$MergedDataSources_Short <- bmassOutput$MergedDataSources[1:1000000,]                                                                                                                                             > system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources_Short, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
     user    system   elapsed
12747.896     5.647 12754.488
> GWASsnps1 <- GWASsnps[1,]
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs2(bmassOutput$MergedDataSources_Short, GWASsnps1, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
Error in GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  incorrect number of subscripts on matrix
In addition: Warning messages:
1: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
2: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
3: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
Timing stopped at: 68.605 0 68.611
> GWASsnps2 <- GWASsnps[1:2,]  
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs2(bmassOutput$MergedDataSources_Short, GWASsnps1, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
   user  system elapsed
105.787   0.000 105.794
Warning messages:
1: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
2: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
3: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
4: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
5: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs2(bmassOutput$MergedDataSources_Short, GWASsnps2, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
   user  system elapsed
212.191   0.003 212.210
Warning messages:
1: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
2: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
3: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
4: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
5: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
6: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
7: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
8: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
9: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
10: In GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) ==  :
  NAs introduced by coercion
.
.
.
> library("lineprof")
> lineprof(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs3(bmassOutput$MergedDataSources_Short, GWASsnps1, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
Reducing depth to 2 (from 6)
Common path:
     time alloc release dups ref                                       src
1   0.006 1.928       0   39  #2 AnnotateDataWithGWASSNPs_Vs2
2  20.633 6.442       0   57  #6 AnnotateDataWithGWASSNPs_Vs2/as.character
3   0.001 0.058       0    0  #6 AnnotateDataWithGWASSNPs_Vs2/as.numeric
4  17.934 1.235       0   47  #6 AnnotateDataWithGWASSNPs_Vs2/as.character
5   0.001 0.523       0    0  #6 AnnotateDataWithGWASSNPs_Vs2/as.numeric
6   0.002 0.012       0   37  #6 AnnotateDataWithGWASSNPs_Vs2
7   2.775 0.475       0   10  #7 AnnotateDataWithGWASSNPs_Vs2/as.character
8   0.001 0.046       0   10  #7 AnnotateDataWithGWASSNPs_Vs2/as.numeric
9  17.923 1.235       0   37  #7 AnnotateDataWithGWASSNPs_Vs2/as.character
10  0.001 0.523       0    0  #7 AnnotateDataWithGWASSNPs_Vs2/as.numeric
11  0.002 0.478       0   37  #7 AnnotateDataWithGWASSNPs_Vs2
12  0.014 0.477       0    0  #9 AnnotateDataWithGWASSNPs_Vs2/==
13  0.001 0.478       0    0  #9 AnnotateDataWithGWASSNPs_Vs2
14  0.014 0.477       0    0 #10 AnnotateDataWithGWASSNPs_Vs2/==
15  0.001 0.000       0    0 #10 AnnotateDataWithGWASSNPs_Vs2
.
.
.
+ #               GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) == GWASsnps1[snpIndex,]$Chr && as.numeric(as.character(MergedDataSource1["BP"])) - BPWindow <= GWASsnps1[snpIndex,]$BP && as.numeric(as.character(MergedDataSource1["BP"])) + BPWindow >= GWASsnps1[snpIndex,]$BP] <- 2
+ #               GWASannot2[as.numeric(as.character(MergedDataSource1["Chr"])) == GWASsnps1[snpIndex,]$Chr && as.numeric(as.character(MergedDataSource1["BP"])) == GWASsnps1[snpIndex,]$BP] <- 1                
+		  GWASannot2[MergedDataSource1["Chr"] == GWASsnps1[snpIndex,]$Chr && MergedDataSource1["BP"] - BPWindow <= GWASsnps1[snpIndex,]$BP && MergedDataSource1["BP"] + BPWindow >= GWASsnps1[snpIndex,]$BP] <- 2
+                 GWASannot2[MergedDataSource1["Chr"] == GWASsnps1[snpIndex,]$Chr && MergedDataSource1["BP"] == GWASsnps1[snpIndex,]$BP] <- 1
.
.
.
> lineprof(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs3(bmassOutput$MergedDataSources_Short, GWASsnps1, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
Reducing depth to 2 (from 6)
Common path:
   time  alloc release dups ref                             src
1 0.002  1.522       0   39  #4 AnnotateDataWithGWASSNPs_Vs2
2 0.278 56.755       0   32  #8 AnnotateDataWithGWASSNPs_Vs2/==
3 0.001  0.478       0    0  #8 AnnotateDataWithGWASSNPs_Vs2
4 0.014  0.477       0    0 #10 AnnotateDataWithGWASSNPs_Vs2/==
5 0.001  0.478       0    0 #10 AnnotateDataWithGWASSNPs_Vs2
6 0.014  0.477       0    0 #11 AnnotateDataWithGWASSNPs_Vs2/==
7 0.001  0.000       0    0 #11 AnnotateDataWithGWASSNPs_Vs2
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs3(bmassOutput$MergedDataSources_Short, GWASsnps1, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
   user  system elapsed
  0.605   0.000   0.605
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs_Vs3(bmassOutput$MergedDataSources_Short, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
   user  system elapsed
 44.119   0.153  44.275
.
.
.
+ #                        if (GWASsnps1[snpIndex,]$BP == as.numeric(as.character(MergedDataSource1["BP"]))) {
+                         if (GWASsnps1[snpIndex,]$BP == MergedDataSource1["BP"]) {
+                                 GWASannot1 <- 1
+                         }
+ #                        else if ((GWASsnps1[snpIndex,]$BP >= as.numeric(as.character(MergedDataSource1["BP"])) - BPWindow) && (GWASsnps1[snpIndex,]$BP <= as.numeric(as.character(MergedDataSource1["BP"])) + BPWindow) && (GWASannot1 != 1)) {
+                         else if ((GWASsnps1[snpIndex,]$BP >= MergedDataSource1["BP"] - BPWindow) && (GWASsnps1[snpIndex,]$BP <= MergedDataSource1["BP"] + BPWindow) && (GWASannot1 != 1)) {
+                                 GWASannot1 <- 2
.
.
.
> system.time(bmassOutput[c("MergedDataSources_Short", "LogFile")] <- AnnotateMergedDataWithGWASSNPs(bmassOutput$MergedDataSources_Short, GWASsnps, GWASsnps_AnnotateWindow, bmassOutput$LogFile)[c("MergedDataSources", "LogFile")])
          ChrBP             Chr              BP             MAF              A1 
  "10_10000135"            "10"     " 10000135"         "0.551"             "g" 
         HDL_A2   HDL_Direction      HDL_pValue           HDL_N      HDL_ZScore 
            "a"             "+" " 4.110000e-01"         "99150"  "8.221351e-01" 
  LDL_Direction      LDL_pValue           LDL_N      LDL_ZScore    TG_Direction 
            "-"    " 9.712e-01"         "94704" "-3.610329e-02"             "-" 
      TG_pValue            TG_N       TG_ZScore    TC_Direction       TC_pValue 
   " 4.344e-01"         "95848" "-7.816845e-01"             "-"    " 8.399e-01" 
           TC_N       TC_ZScore       GWASannot 
       " 99434" "-2.020214e-01"             "0" 
Error in MergedDataSource1["BP"] - BPWindow : 
  non-numeric argument to binary operator
Timing stopped at: 14.4 0.183 14.585 

```







