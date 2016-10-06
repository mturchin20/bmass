# Version News and Updates

Version updates will be tracked and explained here. Major updates & releases will be particularly highlighted.

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

