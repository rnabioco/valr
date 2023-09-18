# valr 0.7.0

* `read_bed` and related functions now automatically calculate the fields. Use of `n_fields` was deprecated.

# valr 0.6.8

* `bed_closest()` now reports all x intervals, even when there are no closest y intervals (e.g. when there is no matching chromosome in y intervals). These intervals are returned populated with `NA` for `.overlap`, `.dist` and y interval locations. 

* Reimplemented `bed_closest()` to use binary search rather than an interval tree search. The closest y interval can be missed with the previous search strategy in high depth interval trees. 

* Fix off by one error when using `max_dist` argument in `bed_cluster()` (#401).

# valr 0.6.7

* Removed `SystemRequirements` from DESCRIPTION to eliminate a NOTE on CRAN.

* `bed_coverage()` now reports intervals from `x` with no matching group in `y` (#395).  

# valr 0.6.6

* Updated [intervalTree](https://github.com/ekg/intervaltree) header to commit f0c4046

* valr now uses [cli](https://cli.r-lib.org/index.html) for more consistent
  errors and messages during interactive use.

* deprecated `genome` argument to `bed_makewindows()` was removed.

# valr 0.6.5

* Handle `max_dist` for first intervals in `bed_cluster()` (#388) 

# valr 0.6.4

* Fixed intron score numbering error in `create_introns` (#377 @sheridar) 

* Fixed bug in handling of list inputs for `bed_intersect()`(#380 @sheridar)   

* Added `read_bigwig` and `read_gtf` functions to import data into valr compatible tibbles (#379)  

* Kent Riemondy is now maintainer.

# valr 0.6.3

* Update to prepare for readr 2.0.0

# valr 0.6.2

## Minor changes

* `RMariaDB` has replaced the deprecated `RMySQL` package as the database backend. 

* valr now imports Rcpp, which should have always been the case,
but was masked by its Import by readr, which recently dropped use of Rcpp.

# valr 0.6.1

## Bug Fixes

* Fixed rchk unprotect error (#365)

# valr 0.6.0

## Major changes

* `trbl_interval()` and `trbl_genome()` custom `tibble` subclasses have been deemed unnecessary and have been removed from the package. 

* coercing `GRanges` to a `valr` compatible data.frame now uses the `gr_to_bed()` function rather than `as.trbl_interal()` methods. 

## Minor changes

* dplyr version < 0.8.0 is no longer supported due to unnecessary code bloat and challenges with handling multiple grouping structures (#359).

* The `sort_by` argument of `bed_random()` has been changed to `sorted`, and will now by default
use `bed_sort()` to sort the output, rather than rely on naming the sorting columns. Sorting can
be suppressed by using `sorted = FALSE`. 

* `bed_sort()` now uses base R sorting with the `radix` method for increased speed. (#353)

* `tbls` processed by `bed_merge()`or `bed_sort()` no longer store either `merged` or `sorted` as attributes, due to these attributes being rarely checked in the codebase and potential sources of unexpected behavior.

## Bug fixes

* Fixed `bed_closest()` to prevent erroneous intervals being reported when adjacent closest intervals are present in the `y` table. (#348)

* Factor columns that are not used for grouping are returned as factors rather than inappropriately being coerced to integer vectors (#360)

# valr 0.5.0

## Major changes 

* Internal `Rcpp` functions have been reorganized to remove all dependencies on `dplyr` C++ functions. 

## Minor changes

* Due to internal refactoring of Rcpp functions, only data.frames containing Numeric, Logical, Integer, Character, and List column types are supported. Columns containing Raw, Complex, or other R classes are not supported and will issue an error. 

* Factors are now disallowed from grouping variables in multiset operations to avoid sort order discrepancies, and compatibility with factor handling in `dplyr` v.0.8.0. Factors will now be internally type-converted to character and a warning is issued.

# valr 0.4.2

## Bug fixes

* Changed the behavior of `as.tbl_interval()` to call `as_tibble()` only on non-tibble input, which prevents groups from being stripped from `tibble()` input (#338).

# valr 0.4.1

* Added new function, `bed_partition()`, which is similar to `bed_merge()` but collapses intervals to elemental intervals rather than the maximal overlapping region. `bed_partition()` also can compute summaries of data from overlapping intervals. See examples in `bed_partition()` and timings in `vignette('benchmarks')` @kriemo.

* Several explicit comparisons to the Bioconductor GenomicRanges library are included for users considering using valr. See examples in `as.tbl_interval()` and timings in `vignette('benchmarks')`.

# valr 0.4.0

## Minor changes

* All relevant tests from bedtool2 were ported into valr. Bugs identified in corner cases by new tests were fixed (#328 @raysinesis)

* `bed_jaccard()` now works with grouped inputs (#216)

* Update dplyr header files to v0.7

* `bed_intersect()` and internal `intersect_impl` were refactored to enable return of non-intersecting intervals.

* The genome argument to `bed_makewindows()` was deprecated and will produce a warning if used. Also error handling was added to check and warn if there are intervals smaller than the requested window size in `makewindows_impl()` (#312 @kriemo)

## Bug fixes

* Fixed off by one error in reported distances from `bed_closest()`. Distances reported now are the same as `bedtools closest` behavior (#311).

* `bed_glyph()` accepts `trbl_intervals` named other than `x` and `y` (#318).

* `bed_makewindows()` now returns the number of windows specified by `num_win` when the input intervals are not evenly divisble into `num_win`, consistent with `bedtools` behavior.

* The output of `findOverlaps()` is now sorted in `subtract_impl()` to prevent reporting intervals that should have been dropped when calling `bed_subtract()` (#316 @kriemo)

# valr 0.3.1

## Enhancements

* A manuscript describing valr has been published in [F1000Research](https://f1000research.com/articles/6-1025/v1).

* New S3 generic `as.tbl_interval()` converts `GenomicRanges::GRanges` objects to `tbl_interval`.

* New `create_tss()` for creating transcription start sites.

* Improve documentation of interval statistics with more complex examples.

## Minor changes

* `bed_sort()` has been de-deprecated to reduce `arrange` calls in library code.

## Bug fixes

* `bed_merge()` now reports start/end columns if spec is provided (#288)

# valr 0.3.0

## Enhancements

* New `create_introns()`, `create_utrs5()` and `create_utrs3()` functions for generating features from BED12 files.

* Speed-ups in `bed_makewindows()` (~50x), `bed_merge()` (~4x), and `bed_flank()` (~4x) (thanks to @kriemo and @sheridar). Thanks to the sponsors of the Biofrontiers Hackathon for the caffeine underlying these improvements.

## Bug fixes

* intervals from `bed_random()` are now sorted properly.

# valr 0.2.0

## Major changes

* Package dplyr v0.5.0 headers with valr to remove dplyr LinkingTo dependency.

* `bed_intersect()` now accepts multiple tbls for intersection (#220 @kriemo).

* new `tbl_interval()` and `tbl_genome()` that wrap tibbles and enforce strict column naming. `trbl_interval()` and `trbl_genome()` are constructors that take `tibble::tribble()` formatting and `is.tbl_interval()` and `is.tbl_genome()` are used to check for valid classes.

## Minor changes

* intervals returned from `bed_random()` are sorted by `chrom` and `start` by default.
  
## Bug fixes

* Merge intervals in `bed_jaccard()` and use numeric values for calculation (fixes #204).

# valr 0.1.2

## Major changes

* Deprecate `bed_sort()` in favor of using `dplyr::arrange()` explicitly (fixes #134).

## Minor changes

* add `src/init.c` that calls `R_registerRoutines` and `R_useDynamicSymbols` to address NOTE in r-devel

* Deprecate `dist` parameter in `bed_closest()` in favor of using user supplied functions (#182 @kriemo)

* Make `.id` values sequential across chroms in `bed_cluster()` output (#171)

* Transfer repository to http://github.com/rnabioco/valr, update links and docs.

* Move shiny app to new repo (http://github.com/rnabioco/valrdata).

* Add Kent Riemondy to LICENSE file.

## Bug fixes

* `bed_merge()` now merges contained intervals (#177)

# valr 0.1.1

## Minor changes

- test / vignette guards for Suggested RMySQL

- fixed memory leak in absdist.cpp

- fixed vignette entry names

# valr 0.1.0

## Major changes

- initial release on CRAN
