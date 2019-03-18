# Version News and Updates

Version updates will be tracked and explained here. Major updates & releases will be particularly highlighted.

## bmass v1.0.0

###### Summary
* Finished prep for CRAN submission
* Setup Travis CI for repo
* Finished leftover unit tests
* Moved vignette setup into two introductory examples
* Edits and suggestions from Peter Carbonetto to help with first main release

###### Since previous version (v0.1.3)
* Finished doing necessary prep and checks for CRAN package submission
  * Finished documenting leftover objects and function arguments for `Roxgyen2` & `devtools::check()`
* Setup and connected Travis CI to git repo
* Finished unit tests in `tests/testthat/test-BasicRFunctionality.R`, `test-ResultsFollowupAnalysisAndPlotting.R`, `test-Stephens2013PLoSONE.globallipids.GLCfuncs.R`, and `test-Stephens2013PLoSONE.test.funcs.R`
* Made `Models` output more user-friendly versions
  * Put phenotype names on top of `Models` and directed users in main vignette to utilize `GetModelPriorMatrix()` for a more interpretable version of `ModelPriors`
* Made a new introductory vignette using a small, simulated dataset, and moved the previous one into a 'more advanced' introductory example (referened as a 'real data' example)
* Merged and adopted edits made by Peter Carbonetto, as well as adopted a few of his suggestions
  * Added proper `releases` on the github repo#
  * Created a `pkgdown` website for the github repo
  * Created the new vignette as described above
  * Updated `README.md`
  * Further cleaned up `roxygen2` comments & examples

###### Notes
* Changed inputs of `GetModelPriorMatrix()` in `ResultsFollowupAnalysisAndPlotting.R` so that `SigmaAlphas` is given a default value to begin with and therefore positioned after `LogFile` 

###### Next steps (if applicable)
* Change output of Matthew's code to use logBFs vs. lbf, and align with overall style of bmass wherever else applicable
* Create second vignette on how to run `bmass()` with some data already created (ie a merged dataset, a phenotype correlation matrix, and priors)
  * Include note that, if provided, `ZScoresCorMatrix` should have the same phenotype order as `DataSources` 
* Double-check new `PassChecksForDataSources.R` additions
* Flesh out `CheckMergedDataSources` as an analogue to `CheckIndividualDataSources` for when a merged dataset is provided.
* Include check that all datasets have same A1 (in `MergeDataGWASAnnotateAndGetMarginalSNPs.R`) & that A1 part of `ExpectedColumnNames` 
* Include some of the extra, downstream functions such as `GetModelPriorMatrix()` in the advanced, introductory vignette


## bmass v0.1.3

###### Summary
* Finished cleaning up all main comments leftover in every `*.R` file
* Added in `Roxygen2` comments and the first vignette, and finished creating unit tests
* Went through `devtools::check()` results and made required edits

###### Since previous version (v0.1.2)
* Cleaned comments & formatting from `bmass.R`, `Stephens2013PLoSONE.globallipids.GLCfuncs.R`, & `Stephens2013PLoSONE.test.funcs.R`
* In `bmass.R` changed `GWASThreshFlag` from `0/1` to `FALSE/TRUE`
  * Mmade related changes to `DetermineAndApplyPriors` in `GetResultsFromMarginalSNPsAndFormat.R`
  * Switched `GWASThreshFlag` to default as `TRUE`
* In `ResultsFollowupAnalysisAndPlotting.R` removed a large number of the functions that were potentially more manuscript-specific and likely to become deprecated (most of which were plotting-related)
* Moved any `library()` calls to proper use of `Imports` in `DESCRIPTION` file

###### Notes
* Dealt with CHECK_0's in `GetResultsFromMarginalSNPsAndFormat.R` (ie `CHECK_0 -- Prob: Throw an error/fail thing here`)
* In `bmass.R` changed `PrintLogStatements` to `PrintProgress` 
* In `bmass.R` removed the loadings of `library(ggplot2)` and `library(reshape2)` 
* In `bmass.R` removed the commented out `ExploreBestModelsUsingPosteriors()`
* In `bmass.R` & main `bmass()` call added `"A1"` to `ExpectedColumnNames`

###### Next steps (if applicable)
* Finish doing necessary prep and checks for CRAN package submission
* Setup and connect Travis CI to git repo
* Change output of Matthew's code to use logBFs vs. lbf, and align with overall style of bmass wherever else applicable
* Fix `Models` and `ModelPriors` outputs to more user-friendly versions
  * Put phenotype names on top of `Models` and possibly condense `ModelPriors` across all the SigmaAlphas
* Create second vignette on how to run `bmass()` with some data already created (ie a merged dataset, a phenotype correlation matrix, and priors)
  * Include note that, if provided, `ZScoresCorMatrix` should have the same phenotype order as `DataSources` 
* Finish unit tests in `test-ResultsFollowupAnalysisAndPlotting.R`, particularly evaluate issues with `GetTopModelsPerSNPViaPosteriors`
* Finish unit tests in `test-Stephens2013PLoSONE.globallipids.GLCfuncs.R` and `test-Stephens2013PLoSONE.test.funcs.R`
* Finish `merge()` unit tests in `tests/testthat/test-BasicRFunctionality.R`
* Flesh out `CheckMergedDataSources` as an analogue to `CheckIndividualDataSources` for when a merged dataset is provided.
* Double-check new `PassChecksForDataSources.R` additions
  * Put in new checks for `Marker` column, input `MergeDataSources` file, etc... as well
* Include check that all datasets have same A1 (in `MergeDataGWASAnnotateAndGetMarginalSNPs.R`) & that A1 part of `ExpectedColumnNames` 
* Check/correct how `write(...stderr())` is functioning in `bmass.R` (since it seems like that might not be occuring like I was expecting it to?) 


## bmass v0.1.2

###### Summary
* Revised some important portions of data-preprocessing
* Began creating and flushing out `ResultsFollowupAnalysisAndPlotting.R` (eg `GetModelPriorMatrix`, `GetTopModelsPerSNPViaPosteriors`, ``, etc...)
* Began cleaning up and adding in proper comments/formatting/etc for eventual, formal, first release (ie v1.0.0)
* Followed up on a number of previous suggested comments/edits from prior work
* Added in `GWASThreshFlag` & `GWASThreshValue` functionality

###### Since previous version (v0.1.1)
* Separated `MergeDataSources()` from `AnnotateMergedDataWithGWASSNPs()` in main `bmass()` call since the two parts take very different times, ie former is longer than the latter.
* Changed `AnnotateDataWithGWASSNPs` to a quicker version where the loop is over the GWAS SNPs and not the merged dataset; taking advantage of in-line conditional statements as well for row selection.
* Included new flag `PrintMergedData` to indicate whether `bmassOutput$MergedDataSources` should be included as part of the output and not made `NULL`. This would be for whether a user wants to save a copy of `MergedDataSources` and reuse it for later runs. 
* Added `bmass()` input variable `ZScoresCorMatrix` that allows users to manually use a specified phenotype correlation matrix.
* Fixed `DetermineAndApplyPriors` to output `PreviousSNPs$logBFs` in the sections when `ProvidedPriors` is not `NULL` and when `UseFlatPriors` is `TRUE`.
* Added in `Marker` column as an expectation for input files in lieu of the `ChrBP` column I constructed. If `ChrBP` is desired to be used as marker names, it can be created by users before submitting the files to `bmass()`. 
* Added in functionality to designate `GWASsnps` not only by input `GWASsnps` file but also by a GWAS p-value threshold (via `GWASThreshFlag` & `GWASThreshValue`)
* Properly designated `PreviousSNPs` sections in `GetResultsFromMarginalSNPsAndFormat.R` (both `DetermineAndApplyPriors` and `FinalizeAndFormatResults`) with `...if (!is.null(GWASsnps)) {...` encapsulation as needed
* Added in `write(..., stderr())` statements to mimic calls to `LogFile1` to be produced as `bmass()` runs, eg giving user a live look into how/where the script is proceeding. These statements are a bit more broad/higher-level than the `LogFile1` calls themselves. Also created `PrintLogStatements` (default `FALSE`) to turn this functionality on/off.
* In `MergeDataSources`, added in section post-merge that checks SNPs still have correct MAF (eg <= .5).
  * Also now checking and dropping SNPs that appear 'fixed' in dataset, eg MAF is equal to 0.
* Moved MAF checks (eg if MAF <= .5, if fixed (0 or 1), etc...) to `CheckIndividualDataSources` section.
  * Added in `CheckDataSourceMAFIsMAF` and `CheckDataSourceMAFFixed` functions.
  * Additionally, now just quitting program if there are fixed SNPs or SNPs that are still > .5; rather have the user fix their datasets and make sure things are correct rather than doing processing mid-package/functions

###### Notes
* Moved `bmass.devlog.vs1.txt` to `bmass.devlog.vs1.md`.
* Changed adding max * length to just add max in `GetSumAcrossSigmaAlphas_withPriors`.
* Moved to data.frame column calls via `$` from atomic vector style calling in `AnnotateDataWithGWASSNPs` after move away from `apply` function.
* Changed `logBF` to `logBFs`.
* Fixed check for `ProvidedPriors` in `CheckIndividualDataSources` to expect 3 ^ (length(DataSources) * length(SigmaAlphas)), not 3 ^ length(DataSources).
* In `bmass.R` changed `if (is.null(PrintMergedData)) {` to `if (PrintMergedData == FALSE) {`, and changed default value of `PrintMergedData` from `NULL` to `FALSE`.
* Changed `posteriorprob` so that `prior` in `apply(prior * 10^lbf,2,normalize)` now specifies the proper matrix dimensions rather than assuming, given proper length, prior will be multiplied along lbf filling column-wise first (which it does by default).
* Added `MergedDataSources` as one of the return values of `ProcessMergedAndAnnotatedDataSources` in `MergeDataGWASAnnotateAndGetMarginalSNPs.R`
* Added `SigmaAlphas` as an input to `GetModelPriorMatrix` in `ResultsFollowupAnalysisAndPlotting.R`
* Corrected `ZScoreHitFlag1` cutoff by including `2*pnorm...`
* Added `if (nrow(SummaryOfTopModels) > 1) {...` to `GetTopModelsPerSNPViaPosteriors` in `ResultsFollowupAnalysisAndPlotting.R` because (apparently) reordering a single-row matrix in the way that I'm doing it messes up its classification as a matrix; after this command, and loss of matrix class, the following `colnames()` command throws an error. This does not occur when there are more than one row and an actualy 're-ordering' occurs (presumably).

###### Next steps (if applicable)
* Flesh out `CheckMergedDataSources` as an analogue to `CheckIndividualDataSources` for when a merged dataset is provided.
* Double-check new `PassChecksForDataSources.R` additions
  * Put in new checks for `Marker` column, input `MergeDataSources` file, etc... as well
* Finish cleaning up comments/edits from offline/pre-dissertation work
  * Start adding in necessary portions for proper R-package structure/eventual CRAN-submission (eg fleshed-out Roxygen2 comments, unit tests, vignettes, necessary additional datafiles for sections such as the unit tests and vignettes, etc...) 
* Include check that all datasets have same A1 & that A1 part of `ExpectedColumnNames` 
* Move towards getting v1.0.0 prepared


## bmass v0.1.1

###### Summary
* First revision and reorganization/restructuring of .R files and internal code
* First steps in setting up unit tests via `testthat`

###### Since previous version (v0.1.0)
* Began moving chunks of code from `PrepareData.R` to better designated .R files
* Changed the few leftover instances of `NewHits` to `NewSNPs`
* Moved `NewSNPs` creation section to after creation of `logBFs` and `PosteriorProbabilities` matrices from MarginalHits; moved `NewHits <- NULL` initialization section of entire code
* Moved from `LogFile1` to `LogFile`
* Changed `ProvidedPriors` to be the first option in the if/elseif/else block that also checks `is.null(GWASsnps)...`, etc...
* Began setting up unit tests via `testthat`, including for basic R functionality, `PassChecksForDataSources.R`, and 
* Introduced main bmass input variable `GWASsnps_AnnotateWindow`
* Moved all instances of `MarginalHits` to `MarginalSNPs`, including in variable names, eg `PruneMarginalHits` and `PruneMarginalHits_bpWindow`
* Split up `GetLogBFsFromData` into `GetLogBFsFromData` and `DetermineAndApplyPriors`
* Moved `GetResultsFromMarginalSNPs.R` to `GetResultsFromMarginalSNPsAndFormat.R`, moved `FinalizeAndFormatResults()` back to `GetResultsFromMarginalSNPsAndFormat.R`, and removed `FormatAndPrepareResults.R`
* Created `ExploreBestModelsUsingPosteriors` to put code chunks dealing with posteriors and related analyses such as `BestClass`
* Moved `FollowupResultsAnalysisFormattingAndPlotting.R` to `ResultsFollowupAnalysisAndPlotting.R`
* Moved `ExploreBestModelsUsingPosteriors` into `ResultsFollowupAnalysisAndPlotting.R`

###### Notes
* Moved `TestData1.txt`, `TestData2.txt`, and `TestSigSNPs.txt` to `bmass_TestData1.txt`, `bmass_TestData2.txt`, and `bmass_TestSigSNPs.txt`
* Currently unclear whether Roxgyen-based .R files describing above test datasets should be in /data or in /R. R CMD check throws warning/error if .R files not in /data, but R won't load datasets via data() properly if .R files not found in /R
* Changed `CheckCharacterFormat` to `CheckCharacterClass`
* Changed underlying code for `CheckDataSourceDirectionColumn`, moving from a for loop to use of `...if (length(...[...!= "+" & ... != "-"]) > 0 )...` 
* Changed `LogFile` to `bmassOutput$LogFile` in main bmass.R file
* Using `bmassOutput` as a way to hold temp variables since can output and assign multiple variables via the list structure and self-made functions. Eg cannot do `c(var1, var2) <- function(input1)`, but can do `list1[c("var1", "var2")] <- function(input1)` assuming function outputs a list with two entries
* Reworked order of input variables for main bmass function
* Included `paste(..., collapse=" ")` at end of `stop()` calls that had `DataSources[!DataSources...]`-type error messages
  * `stop()` seems to already have some automatic `paste()`-type calls? Eg how it can incorporate my call to `paste()` at the end without any additional specifications
* Moved `LogBF` to `logBF` in all instances except for function name `GetLogBFsFromData`
* Moved `if (!is.null(bmassSeedValue)) {...` section to within the `else{...}` portion of the `if (!is.null(ProvidedPriors)) {...` block
* Changed `Priors_Used` to `ModelPriors_Used`
* In `GetSumAcrossSigmaAlphas_withPriors` changed `ModelPriors` to `ModelPriors_Matrix`
* In function call for `GetSumAcrossSigmaAlphas_withPriors` in `FinalizeAndFormatResults`, including `...matrix(...nrow=length(ModelPriors)...)...` to make the number of rows explicit

###### Next steps (if applicable)
* Finish creating unit tests for all .R files and necessary functions
* Finish Roxygen commenting
* Continue revising/reviewing code, including working through remaining CHECK_0 flags
* Go through R CMD CHECK output?


## bmass v0.1.0

###### Summary
* Reworked format of `bmassOutput` to include the `list()` structure of `MaringalSNPs`, `NewSNPs`, and `PreviousSNPs` 

###### Since previous version (v0.0.4)
* Added the following to `bmassOutput`: `NewSNPs`, `NewSNPs$SNPs`, `NewSNPs$logBFs`, `NewSNPs$Posteriors`
  * These are mostly focused on including the GWAS snp results
* Initialized `bmassOutput` variables and associated lists at beginning of function 
* Changed `SNPMarginalpValThreshold` to `SNPMarginalUnivariateThreshold` and `SNPMarginalMultivariateThreshold`
* Moved `PosteriorProbabilities` to `Posteriors`
* Created `MarginalSNPs`, `NewSNPs`, and `PreviousSNPs` setup with `list()` formatting
  * Seeing redundancy among MarginalSNPs, NewSNPs, and Previous SNPs, lead to this suggestion and change
* Moved back to `PreviousSNPs` from `GWAShits`
* Added `BestModel_Posterior` to `PreviousSNPs` mostly has a reminder/placeholder for the time being to expand on it later for both `NewSNPs` and `MarginalSNPs`

###### Notes
* Decided the creation of the first, fully working `bmassOutput` warranted moving to version 0.1.0
  * Partially thinking should have moved to 0.1.0 from 0.0.3 but turns out there were a few useful edits to make to warrant 0.0.4 in the first place
* Next string of updates and 'small patches' will be geared towards reworking and cleaning up the organization of the R files and associated functions
* Then 0.2.0 on might focuse on creating useful/necessary follow-up analysis functions?
* Include `BestModel` as a variable for every SNP? Include posterior as well as the model itself?

###### Next steps (if applicable)
* Same as 0.0.4
* For follow-up functions: Put `ModelMatrix`, `AllBut1Assoc`, and `MarginalPosteriors` into this area?


## bmass v0.0.4

###### Summary
* Created first draft of preliminary, complete `bmassOutput` variable

###### Since previous version (v0.0.3)
* Fleshed out first `is.null(GWASsnps)` & `else` subsections
* Implemented pruning section and `PruneMarginalHits` flag
* Created `logBFs` and `PosteriorProbabilities` output matrices
  * Initialized `GetSumAcrossSigmaAlphas_withPriors()`
* Created `NewSNPs` output variable, in conjunction with `GWASlogBFMinThreshold`
* Added current desired, associated variables to `bmassOutput`

###### Notes
* Current content associated with `bmassOutput`:
  * `MarginalSNPs` 
  * `ModelPriors`
  * `logBFs`
  * `PosteriorProbabilities`
  * `GWASlogBFMinThreshold`
  * `NewSNPs`

###### Next steps (if applicable)
* Next steps include:
  * cleaning up code (eg superfluous comments, print statements, etc...)
  * moving to a return() ending for the code via `bmassOutput`
  * beginning to map out breaking up code into more definable subsections for better organized .R file subsets
  * reviewing, updating, and changing (if/as necessary), `bmassOutput` associated variables
    * any other additions? include both logBFs and PosteriorProbabilities?
  * begin considering downstream helper functions for follow-up analyses/plots of `bmassOutput` results
 

<!---
## bmass v#.#.#

###### Summary

###### Since previous version (v#.#.#)

###### Notes

###### Next steps (if applicable)
-->

