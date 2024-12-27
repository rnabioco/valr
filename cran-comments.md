## New patch version

* This is a new patch version to address UBSAN/ASAN errors from C code in 
a library used by the rtracklayer bioconductor dependency. The rtracklayer 
package has been removed as a dependency to address these errors. 

## Test environment

* win-builder, R devel
* Windows (on Github Actions), R 4.4.2
* macOS (on Github Actions),  R 4.4.2
* ubuntu 22.04.1 (on Github Actions), R devel and 4.4.2
* macOS (local install), R  4.4.1 

## R CMD check results

* on Windows (win-builder, devel)

  Status: 1 NOTE
 * checking CRAN incoming feasibility ... [11s] NOTE
   Maintainer: 'Kent Riemondy <kent.riemondy@gmail.com>'

   Version contains large components (0.8.2.9000)

* on Windows (4.4.2)
  
  Status: OK
  0 errors | 0 warnings | 0 notes
   
* on OS X (R 4.4.1 and 4.4.2)

  Status: OK
  0 errors | 0 warnings | 1 notes
  
  Maintainer: ‘Kent Riemondy <kent.riemondy@gmail.com>’
  
  Version contains large components (0.8.2.9000)
  

* on ubuntu (R 4.4.1)
  
  Status: 1 NOTE
  
  checking installed package size ... NOTE
  installed size is 15.6Mb
  sub-directories of 1Mb or more:
    libs  14.4Mb

* on ubuntu (R devel)  

  Status: OK
  
  checking installed package size ... INFO
  installed size is 15.6Mb
  sub-directories of 1Mb or more:
    libs  14.4Mb
    
## revdepcheck results

We checked 2 reverse dependencies (1 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
