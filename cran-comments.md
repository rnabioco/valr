## New patch version

* This is a new patch version that preemptively fixes an error raised
by the forthcoming readr release.

## Test environment

* win-builder (devel and  R 3.6.0)
* local OS X install, R 3.6.0
* Windows (on Github Actions), 3.6.0
* macOS (on Github Actions), R 3.6.0
* ubuntu 16.04 (on Github Actions), (devel and R 3.6.0)

## R CMD check results

* on rhub (windows, release and devel versions)

  Status: 1 Note
  
  Packages suggested but not available for checking:
  'GenomicRanges' 'IRanges' 'S4Vectors'
  
  GenomicRanges, IRanges, and S4Vectors are available as Bioconductor packages.
  
* on win-builder (release and devel)

  Status: OK
  0 errors | 0 warnings | 1 note
  
  installed size is 18.5Mb
  sub-directories of 1Mb or more:
    libs   17.1Mb
    
  This package uses Rcpp, which creates a large shared library on windows.
  
* on OS X 

  Status: OK
  0 errors | 0 warnings | 1 note
  
  installed size is 18.5Mb
  sub-directories of 1Mb or more:
    libs   17.1Mb
    
  This package uses Rcpp, which creates a large shared library on osx.
  
* on ubuntu

  Status: 1 NOTE
  
  installed size is 18.5Mb
  sub-directories of 1Mb or more:
    libs   17.1Mb

  This package uses Rcpp, which creates a large shared library on linux.
  
## Reverse dependencies

There are no reverse dependencies.
