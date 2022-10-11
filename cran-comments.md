## New patch version

* This is a new patch version that updates a third-party C++ library included in 
the package and adds error reporting via the cli R package.

## Test environment

* win-builder (R-devel and  R 4.2.1)
* local OS X install, R 4.2.0 
* Windows (on Github Actions), 4.2.1 (2022-06-23 ucrt)
* macOS (on Github Actions), R  4.2.1
* ubuntu 20.04.3 (on Github Actions), (devel and R  4.2.1)

## R CMD check results


* on win-builder (release and devel)

  Status: 1 NOTE
  0 errors | 0 warnings | 1 note
  
  * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Kent Riemondy <kent.riemondy@cuanschutz.edu>'

    Version contains large components (0.6.4.9000)
  
    - This package uses Rcpp, which creates a large shared library on windows.
  
* on OS X 

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu (release and devel)

  Status: 1 NOTE
  
  * checking installed package size ... NOTE
    installed size is 19.6Mb
    sub-directories of 1Mb or more:
      libs  18.5Mb

  - This package uses Rcpp, which creates a large shared library on linux.
  
## Reverse dependencies

We checked 1 reverse dependencies (0 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

