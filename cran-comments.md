## New patch version

* This is a new patch version to ensure compatibility with an upcoming release of ggplot 3.6.0. 

## Test environment

* win-builder, R devel
* Windows (on Github Actions), R 4.5.1
* macOS (on Github Actions),  R 4.5.1
* ubuntu 22.04.1 (on Github Actions), R devel and 4.5.1
* macOS (local install), R  4.5.1 

## R CMD check results

* on Windows (win-builder, devel)

  Status: 1 NOTE
* checking CRAN incoming feasibility ... [16s] NOTE
Maintainer: 'Kent Riemondy <kent.riemondy@gmail.com>'

Version contains large components (0.8.3.9000)

* on Windows (4.5.1)
  
  Status: OK
  0 errors | 0 warnings | 0 notes
   
* on OS X (R 4.5.1)

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu (R 4.5.1)
  
  Status: 1 NOTE
  
* checking installed package size ... INFO
  installed size is 15.1Mb
  sub-directories of 1Mb or more:
    libs  13.9Mb

* on ubuntu (R devel)  

  Status: OK
  
  checking installed package size ... INFO
  installed size is 15.1Mb
  sub-directories of 1Mb or more:
    libs  13.9Mb
    
## revdepcheck results

We checked 3 reverse dependencies (2 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 1 packages

Issues with CRAN packages are summarised below.

### Failed to check

* gap (NA)
