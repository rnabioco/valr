## New patch version

* This is a new patch version that fixes a minor bug in the bed_coverage() function, and removes C++11 specifications from the `SystemRequirements` field and Makevar files.   

## Test environment

* win-builder (R-devel)
* local OS X install, (4.2.2) 
* Windows (on Github Actions), 4.2.2 (2022-10-31 ucrt)
* macOS (on Github Actions), R 4.2.2
* ubuntu 22.04.1 (on Github Actions), (devel and R 4.2.2)

## R CMD check results

* on win-builder (devel)

  Status: 1 NOTE
  0 errors | 0 warnings | 1 note
  
  * checking CRAN incoming feasibility ... [18s] NOTE
  Maintainer: 'Kent Riemondy <kent.riemondy@cuanschutz.edu>'

  Version contains large components (0.6.6.9000)

* on windows github actions (R 4.2.2)

  0 errors | 0 warnings | 0 notes
  
* on OS X (R 4.2.2 and devel)

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu (devel)

  Status: 1 NOTE
  
  * checking installed package size ... NOTE
  installed size is 14.2Mb
  sub-directories of 1Mb or more:
    libs  13.0Mb

  - This package uses Rcpp, which creates a large shared library on linux.
  
* on ubuntu (R 4.2.2)

  Status: 1 NOTE
  
  * checking installed package size ... NOTE
  installed size is 14.2Mb
  sub-directories of 1Mb or more:
    libs  13.0Mb

  - This package uses Rcpp, which creates a large shared library on linux.
  

## Reverse dependencies

We checked 1 reverse dependencies (0 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

