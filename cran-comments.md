## New patch version

* This patch release includes bug fixes and performance improvements:
  - `bed_slop()` and `bed_flank()` now preserve input row order
  - Fixed `bed_closest()` to respect custom `suffix` parameter
  - Improved memory efficiency in `bed_intersect()`

* Updated maintainer e-mail address (same maintainer, new generic e-mail)

## Test environment

* win-builder, R devel
* Windows (on Github Actions), R 4.5.2
* macOS (on Github Actions),  R 4.5.2
* ubuntu 22.04.1 (on Github Actions), R devel and 4.5.2
* macOS (local install), R  4.5.2

## R CMD check results

* on Windows (win-builder, devel)

Status: NOTE

  - Non-API call to 'ATTRIB' originates from the cpp11 package dependency

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

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 1 packages

Issues with CRAN packages are summarised below.

### Failed to check

* gap (Fortran linker error unrelated to valr; fails identically with CRAN and dev versions)
