
# valr <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/rnabioco/valr/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/rnabioco/valr/actions/workflows/check-standard.yaml)
[![codecov](https://codecov.io/github/rnabioco/valr/graph/badge.svg)](https://app.codecov.io/github/rnabioco/valr)
[![](https://www.r-pkg.org/badges/version/valr)](https://CRAN.R-project.org/package=valr)
<!-- badges: end -->

valr provides tools to read and manipulate genome intervals and signals,
similar to the [BEDtools](https://bedtools.readthedocs.io/en/latest/)
suite.

## Installation

<div class=".pkgdown-release">

``` r
# Install released version from CRAN
install.packages("valr")
```

</div>

<div class=".pkgdown-devel">

``` r
# Install development version from GitHub
# install.packages("pak")
pak::pak("rnabioco/valr")
```

</div>

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

snps <- read_bed(valr_example("hg19.snps147.chr22.bed.gz"))
genes <- read_bed(valr_example("genes.hg19.chr22.bed.gz"))

# find snps in intergenic regions
intergenic <- bed_subtract(snps, genes)
# find distance from intergenic snps to nearest gene
nearby <- bed_closest(intergenic, genes)

nearby |>
  select(starts_with("name"), .overlap, .dist) |>
  filter(abs(.dist) < 5000)
#> # A tibble: 1,047 × 4
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
#> # ℹ 1,037 more rows
```
