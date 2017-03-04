# valr 0.1.2 (unreleased)

## Major changes

* Deprecate `bed_sort()` in favor of using `dplyr::arrange()` explicitly (fixes #134).

## Minor changes

* Make `.id` values sequential across chroms in `bed_cluster()` output (#171)

* Transfer repository to http:://github.com/rnabioco/valr, update links and docs.

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
