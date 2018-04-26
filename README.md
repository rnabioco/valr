
# valr <img src="man/figures/logo.png" align="right" />

[![Build
Status](https://travis-ci.org/rnabioco/valr.svg?branch=master)](https://travis-ci.org/rnabioco/valr)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/rnabioco/valr?branch=master&svg=true)](https://ci.appveyor.com/project/jayhesselberth/valr)
[![Coverage
Status](https://img.shields.io/codecov/c/github/rnabioco/valr/master.svg)](https://codecov.io/github/rnabioco/valr?branch=master)
[![](https://www.r-pkg.org/badges/version/valr)](https://CRAN.R-project.org/package=valr)
[![](https://cranlogs.r-pkg.org/badges/valr?color=FFD700)](https://www.r-pkg.org/pkg/valr)
[![DOI](https://zenodo.org/badge/49370633.svg)](https://zenodo.org/badge/latestdoi/49370633)

The valr package provides tools to read and manipulate genome intervals
and signals, similar to the
[BEDtools](http://bedtools.readthedocs.io/en/latest/) suite. valr
enables analysis in the R/RStudio environment, leveraging modern R tools
in the [tidyverse](http://tidyverse.org) for a terse, expressive syntax.
Compute-intensive algorithms are implemented in
[Rcpp](http://www.rcpp.org)/C++, and many methods take advantage of the
speed and grouping capability provided by
[dplyr](http://dplyr.tidyverse.org). See `vignette(valr)` for more
details.

## Installation

The latest stable version can be installed from CRAN:

``` r
install.packages('valr')
```

The latest development version can be installed from github:

``` r
# install.packages("devtools")
devtools::install_github('rnabioco/valr')
```

## valr Example

Functions in valr have similar names to their BEDtools counterparts, and
so will be familiar to users coming from the BEDtools suite. Unlike
other tools that wrap BEDtools and write temporary files to disk, valr
tools run natively in memory. Similar to
[pybedtools](https://daler.github.io/pybedtools/#why-pybedtools), valr
has a terse syntax:

``` r
library(valr)
library(dplyr)

snps <- read_bed(valr_example('hg19.snps147.chr22.bed.gz'), n_fields = 6)
genes <- read_bed(valr_example('genes.hg19.chr22.bed.gz'), n_fields = 6)

# find snps in intergenic regions
intergenic <- bed_subtract(snps, genes)
# find distance from intergenic snps to nearest gene
nearby <- bed_closest(intergenic, genes)

nearby %>%
  select(starts_with('name'), .overlap, .dist) %>%
  filter(abs(.dist) < 5000)
#> # A tibble: 1,047 x 4
#>    name.x      name.y   .overlap .dist
#>    <chr>       <chr>       <int> <int>
#>  1 rs530458610 P704P           0  2579
#>  2 rs2261631   P704P           0  -268
#>  3 rs570770556 POTEH           0  -913
#>  4 rs538163832 POTEH           0  -953
#>  5 rs190224195 POTEH           0 -1399
#>  6 rs2379966   DQ571479        0  4750
#>  7 rs142687051 DQ571479        0  3558
#>  8 rs528403095 DQ571479        0  3309
#>  9 rs555126291 DQ571479        0  2745
#> 10 rs5747567   DQ571479        0 -1778
#> # ... with 1,037 more rows
```
