## New patch version

* This is a new patch version to address NOTE on CRAN about Rd link targets, and 
to update the maintainer email address. 

## Test environment

* win-builder, R devel
* Windows (on Github Actions), R 4.3.3
* macOS (on Github Actions),  R 4.3.3
* ubuntu 22.04.1 (on Github Actions), R devel and 4.3.3
* macOS (local install), R 4.3.1

## R CMD check results

* on Windows (win-builder, devel)

  Status: 1 NOTE
 * checking CRAN incoming feasibility ... [11s] NOTE
   Maintainer: 'Kent Riemondy <kent.riemondy@gmail.com>'

   Version contains large components (0.8.1.9000)

  New maintainer:
   Kent Riemondy <kent.riemondy@gmail.com>
  Old maintainer(s):
   Kent Riemondy <kent.riemondy@cuanschutz.edu>

* on Windows (4.4.1)
  
  Status: OK
  0 errors | 0 warnings | 0 notes
   
* on OS X (R 4.4.1)

  Status: OK
  0 errors | 0 warnings | 1 notes
  
  Maintainer: ‘Kent Riemondy <kent.riemondy@gmail.com>’
  
  Version contains large components (0.8.1.9000)
  
  New maintainer:
    Kent Riemondy <kent.riemondy@gmail.com>
  Old maintainer(s):
    Kent Riemondy <kent.riemondy@cuanschutz.edu>

* on ubuntu (R 4.4.1 and R devel)

  Status: 1 NOTE
  * checking installed package size ... NOTE
  installed size is 15.6Mb
  sub-directories of 1Mb or more:
    libs  14.4Mb
    
## revdepcheck results

We checked 2 reverse dependencies (1 from CRAN + 1 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages



