# Version News and Updates

Version updates will be tracked and explained here. Major updates & releases will be particularly highlighted.

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

