




~~~~~~~~~


##20160812

Copying in original code from Matthew's multivariate git repo (https://github.com/stephens999/multivariate) with names that signify their original source. Going to use them as a starting point for editing/cleaning the main function files

```
cp -p /Users/mturchin20/Documents/Work/LabMisc/StephensLab/Multivariate/multivariate/test.funcs.R /Users/mturchin20/Documents/Work/LabMisc/StephensLab/bmass/R/Stephens2013PLoSONE.test.funcs.R
cp -p /Users/mturchin20/Documents/Work/LabMisc/StephensLab/Multivariate/multivariate/globallipids/GLCfuncs.R /Users/mturchin20/Documents/Work/LabMisc/StephensLab/bmass/R/Stephens2013PLoSONE.globallipids.GLCfuncs.R
```


##20161006

Including tags onto current and past commits

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

Practicing/working on R CMD BUILD/CHECK functions

```

```


