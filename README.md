# Targeted Learning for the Sample Average Treatment Effect on Treated Units (SATT)

**2017 Atlantic Causal Inference Competition**

Authors: Jonathan Stiles, Chris J. Kennedy, Caleb Miles, Ivana Malenica, Nima Hejazi, Andre Waschka, Alan Hubbard.

Acknowledgements: We thank Susan Gruber for theoretical inspiration and for sharing the source code from her 2016 entry, and Mark van der Laan for helpful discussion.

## Requirements

* R 3.2 or later, R 3.3+ recommended.
* Java JDK for rJava
* R Packages:
    * CRAN: bartMachine, caret, doMC, earth, ggplot2, glmnet, kernlab, mgcv, nnet, randomForest, RhpcBLASctl, rpart, sandwich, xgboost, xtable
    * Github: ecpolley/SuperLearner, ck37/ck37r

## How to run

Minimal

* Make sure java JDK is installed and R can load rJava & bartMachine packages.
* Run setup.R to install other necessary packages: `make setup`
* Modify targeted_learning.R settings at the top of the file if necessary.
* ./targeted_learning.R inputData outfile1 outfile2

Analysis of 2016 or 2017 data

* Unzip [2017 data](http://faculty.chicagobooth.edu/richard.hahn/pre_data.tar.gz) into `inbound/pre_data/`
* Unzip [2016 data](https://drive.google.com/file/d/0B8TUkApaUlsGekFSblJWa25NM1E/edit) into `inbound/data-2016/`
* Run import-2016.R to import the 2016 data: `make import-2016`
* Run test-2016.R to conduct a single test analysis of 2016: `make test-2016`
* Run analyze-2016.R to analyze all 2016 files using targeted_learning.R: `make analyze-2016`
* Run import-2017.R to import the 2017 data.


## Subdirectory layout

* Archive - old source code that is saved for reference.
* Data - working RData files generated during analysis, not tracked via git.
* Exports - exported files (cvs, tsvs, etc.) that are not tracked via git.
* Inbound - input datasets that are not tracked via git.
* Lib - R source code that defines functions; all .R files are loaded.
* Output - log output files from Savio jobs  etc.
* Scripts - shell (BASH) scripts.
* Simulations - simulation studies.
* Tex - tex output
* Writeup - analysis reports and memos
