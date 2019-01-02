## New minor version

* This is a new major version (see NEWS.md).

* Addresses an incompatibility with the forthcoming dplyr v0.8.0 release.

## Test environments

* Windows Server 2008 R2 SP1, 32/64 bit (on rhub) (devel and R 3.5.1)
* Windows Server 2012 R2 x64 (on appveyor), (devel and R 3.5.1)
* win-builder (devel and  R 3.5.1)
* local OS X install, R 3.5.1
* OS X (on travis-ci), R 3.5.0
* ubuntu 14.04 (on travis-ci), (devel and R 3.5.1)


## R CMD check results

* on appveyor

  Status: OK
  0 errors | 0 warnings | 0 notes
  
* on appveyor (devel version)

  Status: 1 NOTE
  
  Package suggested but not available for checking: 'GenomicRanges'  
  
  GenomicRanges is available as a Bioconductor package.  
  
* on rhub (windows, release and devel versions)

  Status: 1 Note
  
  Packages suggested but not available for checking:
  'GenomicRanges' 'IRanges' 'S4Vectors'
  
  GenomicRanges, IRanges, and S4Vectors are available as Bioconductor packages.
  
* on win-builder (release and devel)

  Status: OK
  0 errors | 0 warnings | 0 notes
  
* on OS X 

  Status: OK
  0 errors | 0 warnings | 0 notes
  
* on ubuntu

  Status: 1 NOTE
  
  installed size is 18.5Mb
  sub-directories of 1Mb or more:
    libs   17.1Mb

  This package uses Rcpp, which creates a large shared library on linux.
  This note is not present on OS X or windows (appveyor or win-builder).
  
## Reverse dependencies

There are no reverse dependencies.
