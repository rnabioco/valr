## New minor version

* I received an email from Dr. Brian Ripley asking for fixes "ASAP" to a memory leak in absdist.cpp and to incorporate checks for the existence of the RMySQL package in the vignettes and tests.

## Test environments

* Windows Server 2012 R2 x64 (on appveyor), R 3.3.2
* win-builder (devel and release)
* local OS X install, R 3.3.2
* OS X (on travis-ci), R 3.3.2
* ubuntu 14.04 (on travis-ci), R 3.3.2

## R CMD check results

* on appveyor

  Status: OK
  0 errors | 0 warnings | 0 notes
 
* on win-builder

  Status: OK
  0 errors | 0 warnings | 0 notes
  
* on OS X 

  Status: OK
  0 errors | 0 warnings | 0 notes
  
* on ubuntu

  Status: 1 NOTE
  
  installed size is 10.7Mb
  sub-directories of 1Mb or more:
    libs   9.7Mb

  This package uses Rcpp, which creates a large shared library on linux.
  This note is not present on OS X or windows (appveyor or win-builder).
  
## Reverse dependencies

There are no reverse dependencies.
