## New minor version (0.10.0)

This is a minor release with new features and two intentional breaking changes
that complete previously announced deprecations:

* The default `min_overlap` value changed from `0` to `1` in `bed_intersect()`,
  `bed_coverage()`, `bed_subtract()`, and `bed_window()`, completing a
  deprecation begun in 0.8.0. Book-ended intervals are now excluded by default,
  matching BEDtools; `min_overlap = 0L` restores the previous behavior.

* The long-deprecated `n_fields` argument of `read_bed()` (deprecated since
  0.6.9) is now defunct.

New features include direct, region-scoped reading of bigWig/bigBed files
(local paths and `http(s)://` URLs) within the interval operations.

## Test environments

* win-builder, R devel
* Windows (GitHub Actions), R release
* macOS (GitHub Actions), R release
* Ubuntu 22.04 (GitHub Actions), R devel and release
* macOS (local install), R 4.6

## R CMD check results

0 errors | 0 warnings | 1 note

* The only NOTE is the expected non-API call to 'ATTRIB', which originates from
  the cpp11 dependency (seen on Windows).

## revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and
dev versions of this package.

 * We saw 0 new problems
 * We failed to check 2 packages

### Failed to check

* gap and tepr fail to install in both the CRAN and dev environments due to their own
  C/Fortran compilation, unrelated to valr.
