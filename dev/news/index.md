# Changelog

## valr (development version)

### Breaking changes

- Added `min_overlap` parameter to
  [`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md),
  [`bed_subtract()`](https://rnabioco.github.io/valr/dev/reference/bed_subtract.md),
  [`bed_coverage()`](https://rnabioco.github.io/valr/dev/reference/bed_coverage.md),
  and
  [`bed_map()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md).
  The default `min_overlap = 1L` aligns with bedtools behavior where
  book-ended (adjacent) intervals are **not** considered overlapping.
  Use `min_overlap = 0L` to preserve the previous valr behavior where
  book-ended intervals were treated as overlapping. Currently, calling
  these functions without an explicit `min_overlap` value will emit a
  deprecation warning and use the legacy behavior (`min_overlap = 0L`).
  In a future release, the default will change to `min_overlap = 1L`, so
  users should update their code to explicitly specify the desired
  behavior.

### Major changes

- Migrated C++ backend from Rcpp to cpp11. This modernizes the codebase
  and removes the Rcpp dependency. All interval operations maintain
  identical behavior.

### Minor changes

- Eliminated all global variable dependencies by replacing bare column
  names with explicit `.data[["column"]]` syntax in data manipulation
  operations and
  [`all_of()`](https://tidyselect.r-lib.org/reference/all_of.html) in
  column selection operations.

- Fixed
  [`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md)
  step size calculation when `step_size` parameter is used. Previously,
  overlapping windows stepped by `win_size - step_size` instead of the
  specified `step_size`
  ([\#438](https://github.com/rnabioco/valr/issues/438)).

- Select methods (`tibble`, `tribble`) are now re-exported from the
  `tibble` package.

- [`read_bigbed()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigbed.html)
  is now re-exported from the `cpp11bigwig` package.

## valr 0.8.4

CRAN release: 2025-06-22

- Update a test for compatibility with forthcoming ggplot2 3.6.0
  ([\#431](https://github.com/rnabioco/valr/issues/431))

## valr 0.8.3

CRAN release: 2025-01-11

- [`read_bigwig()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigwig.html)
  now uses cpp11bigwig on CRAN. The `set_strand` param was removed to be
  more consistent with expected bigWig contents.

- [`read_gtf()`](https://rnabioco.github.io/valr/dev/reference/read_gtf.md)
  was deprecated. The rtracklayer package used for this functionality is
  no longer a dependency of valr due to errors from CRAN AddressSantizer
  checks of the UCSC c-library code vendored in rtracklayer.

- valr now depends on R \>= 4.0.0.

## valr 0.8.2

CRAN release: 2024-08-30

- Address NOTE on CRAN about Rd link targets.

- Change maintainer email address.

## valr 0.8.1

CRAN release: 2024-04-22

- Make vdiffr dependency optional during package testing.

## valr 0.8.0

CRAN release: 2024-04-04

- Added
  [`bed_genomecov()`](https://rnabioco.github.io/valr/dev/reference/bed_genomecov.md)
  to compute interval coverage across a genome.

## valr 0.7.0

CRAN release: 2023-09-18

- `read_bed` and related functions now automatically calculate the
  fields. Use of `n_fields` was deprecated.

## valr 0.6.8

CRAN release: 2023-05-16

- [`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md)
  now reports all x intervals, even when there are no closest y
  intervals (e.g. when there is no matching chromosome in y intervals).
  These intervals are returned populated with `NA` for `.overlap`,
  `.dist` and y interval locations.

- Reimplemented
  [`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md)
  to use binary search rather than an interval tree search. The closest
  y interval can be missed with the previous search strategy in high
  depth interval trees.

- Fix off by one error when using `max_dist` argument in
  [`bed_cluster()`](https://rnabioco.github.io/valr/dev/reference/bed_cluster.md)
  ([\#401](https://github.com/rnabioco/valr/issues/401)).

## valr 0.6.7

CRAN release: 2023-02-18

- Removed `SystemRequirements` from DESCRIPTION to eliminate a NOTE on
  CRAN.

- [`bed_coverage()`](https://rnabioco.github.io/valr/dev/reference/bed_coverage.md)
  now reports intervals from `x` with no matching group in `y`
  ([\#395](https://github.com/rnabioco/valr/issues/395)).

## valr 0.6.6

CRAN release: 2022-10-11

- Updated [intervalTree](https://github.com/ekg/intervaltree) header to
  commit f0c4046

- valr now uses [cli](https://cli.r-lib.org/index.html) for more
  consistent errors and messages during interactive use.

- deprecated `genome` argument to
  [`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md)
  was removed.

## valr 0.6.5

CRAN release: 2022-08-19

- Handle `max_dist` for first intervals in
  [`bed_cluster()`](https://rnabioco.github.io/valr/dev/reference/bed_cluster.md)
  ([\#388](https://github.com/rnabioco/valr/issues/388))

## valr 0.6.4

CRAN release: 2021-12-08

- Fixed intron score numbering error in `create_introns`
  ([\#377](https://github.com/rnabioco/valr/issues/377)
  [@sheridar](https://github.com/sheridar))

- Fixed bug in handling of list inputs for
  [`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md)([\#380](https://github.com/rnabioco/valr/issues/380)
  [@sheridar](https://github.com/sheridar))

- Added `read_bigwig` and `read_gtf` functions to import data into valr
  compatible tibbles
  ([\#379](https://github.com/rnabioco/valr/issues/379))

- Kent Riemondy is now maintainer.

## valr 0.6.3

CRAN release: 2021-05-15

- Update to prepare for readr 2.0.0

## valr 0.6.2

CRAN release: 2020-10-07

### Minor changes

- `RMariaDB` has replaced the deprecated `RMySQL` package as the
  database backend.

- valr now imports Rcpp, which should have always been the case, but was
  masked by its Import by readr, which recently dropped use of Rcpp.

## valr 0.6.1

CRAN release: 2020-05-08

### Bug Fixes

- Fixed rchk unprotect error
  ([\#365](https://github.com/rnabioco/valr/issues/365))

## valr 0.6.0

CRAN release: 2020-05-04

### Major changes

- `trbl_interval()` and `trbl_genome()` custom `tibble` subclasses have
  been deemed unnecessary and have been removed from the package.

- coercing `GRanges` to a `valr` compatible data.frame now uses the
  [`gr_to_bed()`](https://rnabioco.github.io/valr/dev/reference/gr_to_bed.md)
  function rather than `as.trbl_interal()` methods.

### Minor changes

- dplyr version \< 0.8.0 is no longer supported due to unnecessary code
  bloat and challenges with handling multiple grouping structures
  ([\#359](https://github.com/rnabioco/valr/issues/359)).

- The `sort_by` argument of
  [`bed_random()`](https://rnabioco.github.io/valr/dev/reference/bed_random.md)
  has been changed to `sorted`, and will now by default use
  [`bed_sort()`](https://rnabioco.github.io/valr/dev/reference/bed_sort.md)
  to sort the output, rather than rely on naming the sorting columns.
  Sorting can be suppressed by using `sorted = FALSE`.

- [`bed_sort()`](https://rnabioco.github.io/valr/dev/reference/bed_sort.md)
  now uses base R sorting with the `radix` method for increased speed.
  ([\#353](https://github.com/rnabioco/valr/issues/353))

- `tbls` processed by
  [`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md)or
  [`bed_sort()`](https://rnabioco.github.io/valr/dev/reference/bed_sort.md)
  no longer store either `merged` or `sorted` as attributes, due to
  these attributes being rarely checked in the codebase and potential
  sources of unexpected behavior.

### Bug fixes

- Fixed
  [`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md)
  to prevent erroneous intervals being reported when adjacent closest
  intervals are present in the `y` table.
  ([\#348](https://github.com/rnabioco/valr/issues/348))

- Factor columns that are not used for grouping are returned as factors
  rather than inappropriately being coerced to integer vectors
  ([\#360](https://github.com/rnabioco/valr/issues/360))

## valr 0.5.0

CRAN release: 2019-01-03

### Major changes

- Internal `Rcpp` functions have been reorganized to remove all
  dependencies on `dplyr` C++ functions.

### Minor changes

- Due to internal refactoring of Rcpp functions, only data.frames
  containing Numeric, Logical, Integer, Character, and List column types
  are supported. Columns containing Raw, Complex, or other R classes are
  not supported and will issue an error.

- Factors are now disallowed from grouping variables in multiset
  operations to avoid sort order discrepancies, and compatibility with
  factor handling in `dplyr` v.0.8.0. Factors will now be internally
  type-converted to character and a warning is issued.

## valr 0.4.2

CRAN release: 2018-11-17

### Bug fixes

- Changed the behavior of `as.tbl_interval()` to call
  [`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  only on non-tibble input, which prevents groups from being stripped
  from [`tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  input ([\#338](https://github.com/rnabioco/valr/issues/338)).

## valr 0.4.1

CRAN release: 2018-06-08

- Added new function,
  [`bed_partition()`](https://rnabioco.github.io/valr/dev/reference/bed_partition.md),
  which is similar to
  [`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md)
  but collapses intervals to elemental intervals rather than the maximal
  overlapping region.
  [`bed_partition()`](https://rnabioco.github.io/valr/dev/reference/bed_partition.md)
  also can compute summaries of data from overlapping intervals. See
  examples in
  [`bed_partition()`](https://rnabioco.github.io/valr/dev/reference/bed_partition.md)
  and timings in `vignette('benchmarks')`
  [@kriemo](https://github.com/kriemo).

- Several explicit comparisons to the Bioconductor GenomicRanges library
  are included for users considering using valr. See examples in
  `as.tbl_interval()` and timings in `vignette('benchmarks')`.

## valr 0.4.0

CRAN release: 2018-01-25

### Minor changes

- All relevant tests from bedtool2 were ported into valr. Bugs
  identified in corner cases by new tests were fixed
  ([\#328](https://github.com/rnabioco/valr/issues/328)
  [@raysinesis](https://github.com/raysinesis))

- [`bed_jaccard()`](https://rnabioco.github.io/valr/dev/reference/bed_jaccard.md)
  now works with grouped inputs
  ([\#216](https://github.com/rnabioco/valr/issues/216))

- Update dplyr header files to v0.7

- [`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md)
  and internal `intersect_impl` were refactored to enable return of
  non-intersecting intervals.

- The genome argument to
  [`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md)
  was deprecated and will produce a warning if used. Also error handling
  was added to check and warn if there are intervals smaller than the
  requested window size in `makewindows_impl()`
  ([\#312](https://github.com/rnabioco/valr/issues/312)
  [@kriemo](https://github.com/kriemo))

### Bug fixes

- Fixed off by one error in reported distances from
  [`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md).
  Distances reported now are the same as `bedtools closest` behavior
  ([\#311](https://github.com/rnabioco/valr/issues/311)).

- [`bed_glyph()`](https://rnabioco.github.io/valr/dev/reference/bed_glyph.md)
  accepts `trbl_intervals` named other than `x` and `y`
  ([\#318](https://github.com/rnabioco/valr/issues/318)).

- [`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md)
  now returns the number of windows specified by `num_win` when the
  input intervals are not evenly divisble into `num_win`, consistent
  with `bedtools` behavior.

- The output of `findOverlaps()` is now sorted in `subtract_impl()` to
  prevent reporting intervals that should have been dropped when calling
  [`bed_subtract()`](https://rnabioco.github.io/valr/dev/reference/bed_subtract.md)
  ([\#316](https://github.com/rnabioco/valr/issues/316)
  [@kriemo](https://github.com/kriemo))

## valr 0.3.1

CRAN release: 2017-07-22

### Enhancements

- A manuscript describing valr has been published in
  [F1000Research](https://f1000research.com/articles/6-1025/v1).

- New S3 generic `as.tbl_interval()` converts
  [`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  objects to `tbl_interval`.

- New
  [`create_tss()`](https://rnabioco.github.io/valr/dev/reference/create_tss.md)
  for creating transcription start sites.

- Improve documentation of interval statistics with more complex
  examples.

### Minor changes

- [`bed_sort()`](https://rnabioco.github.io/valr/dev/reference/bed_sort.md)
  has been de-deprecated to reduce `arrange` calls in library code.

### Bug fixes

- [`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md)
  now reports start/end columns if spec is provided
  ([\#288](https://github.com/rnabioco/valr/issues/288))

## valr 0.3.0

CRAN release: 2017-06-15

### Enhancements

- New
  [`create_introns()`](https://rnabioco.github.io/valr/dev/reference/create_introns.md),
  [`create_utrs5()`](https://rnabioco.github.io/valr/dev/reference/create_utrs5.md)
  and
  [`create_utrs3()`](https://rnabioco.github.io/valr/dev/reference/create_utrs3.md)
  functions for generating features from BED12 files.

- Speed-ups in
  [`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md)
  (~50x),
  [`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md)
  (~4x), and
  [`bed_flank()`](https://rnabioco.github.io/valr/dev/reference/bed_flank.md)
  (~4x) (thanks to [@kriemo](https://github.com/kriemo) and
  [@sheridar](https://github.com/sheridar)). Thanks to the sponsors of
  the Biofrontiers Hackathon for the caffeine underlying these
  improvements.

### Bug fixes

- intervals from
  [`bed_random()`](https://rnabioco.github.io/valr/dev/reference/bed_random.md)
  are now sorted properly.

## valr 0.2.0

CRAN release: 2017-05-05

### Major changes

- Package dplyr v0.5.0 headers with valr to remove dplyr LinkingTo
  dependency.

- [`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md)
  now accepts multiple tbls for intersection
  ([\#220](https://github.com/rnabioco/valr/issues/220)
  [@kriemo](https://github.com/kriemo)).

- new `tbl_interval()` and `tbl_genome()` that wrap tibbles and enforce
  strict column naming. `trbl_interval()` and `trbl_genome()` are
  constructors that take
  [`tibble::tribble()`](https://tibble.tidyverse.org/reference/tribble.html)
  formatting and `is.tbl_interval()` and `is.tbl_genome()` are used to
  check for valid classes.

### Minor changes

- intervals returned from
  [`bed_random()`](https://rnabioco.github.io/valr/dev/reference/bed_random.md)
  are sorted by `chrom` and `start` by default.

### Bug fixes

- Merge intervals in
  [`bed_jaccard()`](https://rnabioco.github.io/valr/dev/reference/bed_jaccard.md)
  and use numeric values for calculation (fixes
  [\#204](https://github.com/rnabioco/valr/issues/204)).

## valr 0.1.2

CRAN release: 2017-03-16

### Major changes

- Deprecate
  [`bed_sort()`](https://rnabioco.github.io/valr/dev/reference/bed_sort.md)
  in favor of using
  [`dplyr::arrange()`](https://dplyr.tidyverse.org/reference/arrange.html)
  explicitly (fixes
  [\#134](https://github.com/rnabioco/valr/issues/134)).

### Minor changes

- add `src/init.c` that calls `R_registerRoutines` and
  `R_useDynamicSymbols` to address NOTE in r-devel

- Deprecate `dist` parameter in
  [`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md)
  in favor of using user supplied functions
  ([\#182](https://github.com/rnabioco/valr/issues/182)
  [@kriemo](https://github.com/kriemo))

- Make `.id` values sequential across chroms in
  [`bed_cluster()`](https://rnabioco.github.io/valr/dev/reference/bed_cluster.md)
  output ([\#171](https://github.com/rnabioco/valr/issues/171))

- Transfer repository to <http://github.com/rnabioco/valr>, update links
  and docs.

- Move shiny app to new repo (<http://github.com/rnabioco/valrdata>).

- Add Kent Riemondy to LICENSE file.

### Bug fixes

- [`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md)
  now merges contained intervals
  ([\#177](https://github.com/rnabioco/valr/issues/177))

## valr 0.1.1

CRAN release: 2016-12-01

### Minor changes

- test / vignette guards for Suggested RMySQL

- fixed memory leak in absdist.cpp

- fixed vignette entry names

## valr 0.1.0

CRAN release: 2016-11-21

### Major changes

- initial release on CRAN
