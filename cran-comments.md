## Resubmission

This is a resubmission. In this version I have:

* corrected a non-canonical CRAN URL.

## New minor version

* Maintenance release to fix some small errors (see NEWS.md)

## Test environments

* Windows Server 2012 R2 x64 (on appveyor), R 3.3.3
* win-builder (devel and release)
* local OS X install, R 3.3.3
* OS X (on travis-ci), R 3.3.3
* ubuntu 14.04 (on travis-ci), R 3.3.3

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
