# valr

[![Build Status](https://travis-ci.org/jayhesselberth/valr.svg?branch=master)](https://travis-ci.org/jayhesselberth/valr)
[![Coverage Status](https://img.shields.io/codecov/c/github/jayhesselberth/valr/master.svg)](https://codecov.io/github/jayhesselberth/valr?branch=master)

`valr` provides tools for fast, exploratory analysis of large biological data sets **within R/RStudio**. `valr` was developed from a desire to teach the powerful concepts of genome interval arithmetic, but without the cumbersome back-and-forth between command-line and exploratory analysis tools. 

Key parts of `valr` are implemented in [`Rcpp`][3] for speed. Moreover, `valr` integrates with [`dplyr`][2] and the `magrittr` pipe operator (`%>%`) for an expressive syntax and can be used in reproducible reports written in [`RMarkdown`][10].

See the [`valr` documentation](http://jayhesselberth.github.io/valr) for a complete API reference and examples.

## Installation

`valr` can be installed from github:

```R
devtools::install_github('jayhesselberth/valr')
```
 
__Note__ that `valr` requires a full-featured C++11 compiler (`gcc>=6.0` or `clang++`)

## API

Function names are similar to their their [BEDtools][1] counterparts, with some additions.

### Reading data

* BED and related files are read with `read_bed()`, `read_bed12()`, `read_bedgraph()`, `read_narrowpeak()` and `read_broadpeak()`.
  
* Genome files containing chromosome name and size information are loaded with `read_genome()`.
  
* VCF files are loaded with `read_vcf()`.

### Transforming single interval sets

* Intervals are ordered with `bed_sort()`.

* Interval coordinates are adjusted with `bed_slop()` and `bed_shift()`, and new flanking intervals are created with `bed_flank()`.

* Nearby intervals are combined with `bed_merge()` and identified (but not merged) with `bed_cluster()`.  

* Intervals not covered by a query are created with `bed_complement()`.

### Comparing multiple interval sets

* Find overlaps between two sets of intervals with `bed_intersect()`.

* Apply functions to selected columns for overlapping intervals with `bed_map()`.

* Remove intervals based on overlaps between two files with `bed_subtract()`.

* Find overlapping intervals within a window with `bed_window()`.

* Find the closest intervals independent of overlaps with `bed_closest()`.

### Randomizing intervals

* Generate random intervals from an input genome with `bed_random()`.

* Shuffle the coordinates of input intervals with `bed_shuffle()`.

* Random sampling of input intervals is done with the `sample_` function family in `dplyr`.

### Interval statistics

* Quality overlaps between two sets of intervals with `bed_fisher()`.

* Quantify relative and absolute distances between sets of intervals with `bed_reldist()` and `bed_absdist()`.

* Quantify extent of overlap between two sets of intervals with `bed_jaccard()`.

## Related work

* Command-line tools [BEDtools][1] and [bedops][5].

* The Python library [pybedtools][4] wraps BEDtools.

* The R packages [GenomicRanges][6], [bedr][7], [IRanges][8] and [GenometriCorr][9] provide similar capability with a different philosophy.

[1]: http://bedtools.readthedocs.org/en/latest/
[2]: https://github.com/hadley/dplyr
[3]: http://www.rcpp.org/
[4]: https://pythonhosted.org/pybedtools/
[5]: http://bedops.readthedocs.org/en/latest/index.html
[6]: https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
[7]: https://cran.r-project.org/web/packages/bedr/index.html
[8]: https://bioconductor.org/packages/release/bioc/html/IRanges.html
[9]: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529
[10]: http://rmarkdown.rstudio.com/
