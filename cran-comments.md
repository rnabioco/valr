## New patch version

* This is a new patch version, with a re-implementation of the C++ internals for `bed_closest()`, and a minor bug fix in `bed_cluster()`.

## Test environment

* win-builder, R devel
* macOS (local install), R 4.3.0
* Windows (on Github Actions), R 4.3.0
* macOS (on Github Actions),  R 4.3.0
* ubuntu 22.04.1 (on Github Actions), R devel and 4.3.0

## R CMD check results

* on Windows (win-builder, devel)

  Status: 1 NOTE
  * checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Kent Riemondy <kent.riemondy@cuanschutz.edu>'

  Version contains large components (0.6.7.9000)

* on Windows (4.3.0)
  
  Status: OK
  0 errors | 0 warnings | 0 notes
   
* on OS X (R 4.3.0)

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu (devel)

  Status: 1 NOTE
  * checking installed package size ... NOTE
  installed size is 14.5Mb
  sub-directories of 1Mb or more:
    libs  13.3Mb
    
* on ubuntu (R 4.3.0)

  Status: 1 NOTE
  * checking installed package size ... NOTE
  installed size is 14.5Mb
  sub-directories of 1Mb or more:
    libs  13.3Mb

## Reverse dependencies

We checked 1 reverse dependencies (0 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

