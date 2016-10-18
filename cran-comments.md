## Test environments

* local OS X install, R 3.3.1
* OS X (on travis-ci), R 3.3.1
* ubuntu 14.04 (on travis-ci), R 3.3.1
* win-builder (devel and release)
* Windows Server 2012 R2 x64 (on appveyor), R 3.3.1

## R CMD check results

* on win-builder and appveyor

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
  This note is not present on OS X or windows.
  
* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
