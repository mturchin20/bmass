## Resubmission
This is a resubmission. In this version I have:

* Edited the Description field in DESCRIPTION to
  match proper, requested formatting and to 
  include M. Stephens reference.

* In the comments from the initial submission it 
  was requested to "provide small executable examples 
  in all of your Rd-files" and "They examples are 
  wrapped in \dontrun{}, hence nothing gets tested."
  However, please note that in each Rd-file where 
  appropriate there is example code outside of the 
  \dontrun{} chunks, and that commands such as 
  'example("bmass")' run properly. Thank you.

## Test environments
* local OS X install, R 3.5.1
* ubuntu 14.04.5 (on travis-ci), R 3.5.2
* win-builder (release & devel), R 3.5.3 & R unstable 
  (2019-05-14 r76503)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE (win-builder):

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Michael Turchin <mturchin20@uchicago.edu>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Phenotypes (22:52)
  Unassociated (17:5)
  phenotypes (15:22)


These words are spelled correctly.
