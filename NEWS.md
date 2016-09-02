# Version News and Updates

Version updates will be tracked and explained here. Major updates & releases will be particularly highlighted.

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
 

