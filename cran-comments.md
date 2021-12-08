## New patch version

* This is a new patch version that provides two new data import functions, fixes for two errors, and names Kent Riemondy as the new maintainer of the package. 

## Test environment

* win-builder (R-devel and  R 4.1.2)
* local OS X install, R 4.1.2 
* Windows (on Github Actions), R 4.1.2
* macOS (on Github Actions), R 4.1.2
* ubuntu 20.04.3 (on Github Actions), (devel and R 4.1.2)

## R CMD check results


* on win-builder (release and devel)

  Status: 1 NOTE
  0 errors | 0 warnings | 1 note
  
  * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Kent Riemondy <kent.riemondy@cuanschutz.edu>'

    Version contains large components (0.6.3.9000)
    
    New maintainer:
      Kent Riemondy <kent.riemondy@cuanschutz.edu>
    
    Old maintainer(s):
      Jay Hesselberth <jay.hesselberth@gmail.com>
  
    - This package uses Rcpp, which creates a large shared library on windows.
  
* on OS X 

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu

  Status: 1 NOTE
  
  * checking installed package size ... NOTE
    installed size is 19.8Mb
    sub-directories of 1Mb or more:
      libs  18.7Mb

  - This package uses Rcpp, which creates a large shared library on linux.
  
## Reverse dependencies

We checked 1 reverse dependencies (0 from CRAN + 1 from Bioconductor (RLSeq v1.0.0)), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

