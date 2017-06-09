# valr 0.2.0.9000

## Enhancements

* New `create_introns()`, `create_utrs5()` and `create_utrs3()` functions for generating features from BED12 files.

* Speed-ups in `bed_makewindows()` (~50x), `bed_merge()` (~4x), and `bed_flank()` (~4x) (thanks to @kriemo and @sheridar). Thanks to the sponsors of the [Biofrontiers Hackathon](https://hackathon.colorado.edu/) for the caffeine underlying these improvements.

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
