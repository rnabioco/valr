## New minor version

* This is a new minor version (see NEWS.md).

* This version incorporates changes to be compatible with the new version of dplyr (v0.7.0).

## Test environments

* Windows Server 2012 R2 x64 (on appveyor), R 3.4.0
* win-builder (devel and release)
* local OS X install, R 3.4.0
* OS X (on travis-ci), R 3.4.0
* ubuntu 14.04 (on travis-ci), R 3.4.0

## R CMD check results

* on appveyor

  Status: OK
  0 errors | 0 warnings | 0 notes
 
* on win-builder

  Status: 1 WARNING
  0 errors | 1 warnings | 0 notes
  
  This WARNING is caused by an issue with pandoc, which I understand CRAN is already addressing.
  
  * checking top-level files ... WARNING
    Conversion of 'README.md' failed:
    pandoc.exe: Could not fetch https://img.shields.io/codecov/c/github/rnabioco/valr/master.svg
    TlsExceptionHostPort (HandshakeFailed Error_EOF) "img.shields.io" 443
  
* on OS X 

  Status: OK
  0 errors | 0 warnings | 0 notes
  
  The build on OS X fails on travis-ci, because the dplyr 0.7.0 binary is not available due to an ERROR in the CRAN binary build. This will be fixed in forthcoming patch release of dplyr (https://github.com/tidyverse/dplyr/commit/0b40356fe22bc33fc9ba38bac22aa898c1c40480)
  
* on ubuntu

  Status: 1 NOTE
  
  installed size is 10.7Mb
  sub-directories of 1Mb or more:
    libs   9.7Mb

  This package uses Rcpp, which creates a large shared library on linux.
  This note is not present on OS X or windows (appveyor or win-builder).
  
## Reverse dependencies

There are no reverse dependencies.
