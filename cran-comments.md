## New minor version

* This is a new minor version that updates the C++ backend from Rcpp to cpp11,
  and updates the maintainer.

## Test environment

* win-builder, R devel
* Windows (on Github Actions), R 4.5.2
* macOS (on Github Actions),  R 4.5.2
* ubuntu 22.04.1 (on Github Actions), R devel and 4.5.2
* macOS (local install), R  4.5.2

## R CMD check results

* on Windows (win-builder, devel)

Status: NOTE

* on Windows (4.5.1)

  Status: OK

* on macos-latest (R 4.5.2)

  Status: OK
  0 errors | 0 warnings | 0 notes

* on ubuntu (R 4.5.2)

  Status: OK

* checking installed package size ... INFO
  installed size is 10.4Mb
  sub-directories of 1Mb or more:
    libs   9.2Mb

* on ubuntu (R devel)

  Status: OK

* checking installed package size ... INFO
  installed size is 10.4Mb
  sub-directories of 1Mb or more:
    libs   9.2Mb

## revdepcheck results

We checked 3 reverse dependencies (2 from CRAN), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 1 new problems
 * We failed to check 1 packages

Issues with CRAN packages are summarised below.

### New problems

* tepr (seems like a false positive, does not build)

### Failed to check

* gap (NA)
