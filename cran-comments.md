## Resubmission
 
This is a resubmission. In this fourth version I:

* excluded README.md from the build to eliminate a 404 error.

In the third (previous) submssion, I:

* updated the DESCRIPTION by singly-quoting and describing the use of 'BEDtools'.

In the second (previous) submission, I:

* Updated to the CRAN MIT LICENSE template.

* Updated the Authors@R section using `person()` calls

* Fixed non-canonical CRAN URLs

* Fixed compilation WARNINGS on CRAN win-builder. These warnings were caused by the use of GNU-specific extensions in Rcpp 0.12.7 (RcppCore/Rcpp#537). The new release of Rcpp (0.12.8) has been updated to eliminate these warnings.

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

  Status: 1 NOTE
  
  The NOTE indicates this is a New submission.
 
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
  
* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
