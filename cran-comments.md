## Resubmission
 
This is a resubmission. In this version I have:

* Updated to the CRAN MIT LICENSE template.

* Updated the Authors@R section using `person()` calls

* Fixed non-canonical CRAN URLs

## Test environments

* Windows Server 2012 R2 x64 (on appveyor), R 3.3.1
* win-builder (devel and release)
* local OS X install, R 3.3.1
* OS X (on travis-ci), R 3.3.1
* ubuntu 14.04 (on travis-ci), R 3.3.1

## R CMD check results

* on appveyor

  Status: OK
  0 errors | 0 warnings | 0 notes
 
* on win-builder

  Status: 1 WARNING, 1 NOTE
  
  The NOTE indicates this is a New submission.
  
  The WARNING is as follows:
  
  Found the following significant warnings:
  d:\RCompile\CRANpkg\lib\3.4/dplyr/include/tools/all_na.h:15:40:
    warning: ISO C99 requires rest arguments to be used
    
  This WARNING is in the dplyr dependency, and is not present in the appveyor / windows build.
 
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
