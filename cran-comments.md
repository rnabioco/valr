## New minor version

* This is a new minor version with limited updates to package functionality.

## Test environment

* macOS (local install), R 4.3.1
* win-builder, R devel
* Windows (on Github Actions), R 4.3.1
* macOS (on Github Actions),  R 4.3.1
* ubuntu 22.04.1 (on Github Actions), R devel and 4.3.1

## R CMD check results

* on Windows (win-builder, devel)

  Status: 1 NOTE
  * checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Kent Riemondy <kent.riemondy@cuanschutz.edu>'

  Version contains large components (0.6.8.9000)

* on Windows (4.3.1)
  
  Status: OK
  0 errors | 0 warnings | 0 notes
   
* on OS X (R 4.3.1)

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu (R 4.3.1)

  Status: 1 NOTE
  * checking installed package size ... NOTE
  installed size is 14.5Mb
  sub-directories of 1Mb or more:
    libs  13.3Mb
    
## Reverse dependencies

We checked 1 reverse dependencies (0 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
