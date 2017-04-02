# valr 0.1.2.9000

## Major changes

* new `tbl_interval()` and `tbl_sizes()` table types that wrap tibbles and enforce strict column naming. Input tables can be and are coerced with `as_tbl_ivl()` and `as_tbl_szs()` in each method.

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
