## Test environments
* local OS X install, R 3.5.3
* ubuntu 14.04.5 (on travis-ci), R 3.5.2
* win-builder (release), R 3.5.3

## R CMD check results
There were no ERRORs, WARNINGs, or NOTES. 

## Downstream dependencies
I have also run R CMD check on downstream dependencies of httr 
(https://github.com/wch/checkresults/blob/master/httr/r-release). 
All packages that I could install passed except:

* Ecoengine: this appears to be a failure related to config on 
  that machine. I couldn't reproduce it locally, and it doesn't 
  with httr 0.4).
