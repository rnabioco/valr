# Compute coverage of intervals.

Compute coverage of intervals.

## Usage

``` r
bed_coverage(x, y, ..., min_overlap = NULL)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- ...:

  extra arguments (not used)

- min_overlap:

  minimum overlap in base pairs required for coverage. Set to `1` to
  exclude book-ended intervals (matching bedtools behavior), or `0` to
  include them (legacy valr behavior). The default will change from `0`
  to `1` in a future version.

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) with
the following additional columns:

- `.ints` number of `x` intersections

- `.cov` per-base coverage of `x` intervals

- `.len` total length of `y` intervals covered by `x` intervals

- `.frac` `.len` scaled by the number of `y` intervals

## Details

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/dev/reference/flip_strands.md).

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html>

Other multiple set operations:
[`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md),
[`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md),
[`bed_map()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md),
[`bed_subtract()`](https://rnabioco.github.io/valr/dev/reference/bed_subtract.md),
[`bed_window()`](https://rnabioco.github.io/valr/dev/reference/bed_window.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 100,    500,  "+",
  "chr2", 200,    400,  "+",
  "chr2", 300,    500,  "-",
  "chr2", 800,    900,  "-"
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value, ~strand,
  "chr1", 150,    400,  100,    "+",
  "chr1", 500,    550,  100,    "+",
  "chr2", 230,    430,  200,    "-",
  "chr2", 350,    430,  300,    "-"
)

bed_coverage(x, y)
#> Warning: The `min_overlap` argument of `bed_coverage()` is deprecated as of valr 0.8.0.
#> ℹ The default will change from 0 (book-ended intervals overlap) to 1 (strict
#>   overlap) in a future version.
#> ℹ Set `min_overlap = 0L` to keep the legacy behavior, or `min_overlap = 1L` for
#>   bedtools-compatible behavior.
#> # A tibble: 4 × 8
#>   chrom start   end strand .ints  .cov  .len .frac
#>   <chr> <dbl> <dbl> <chr>  <int> <int> <dbl> <dbl>
#> 1 chr1    100   500 +          2   250   400 0.625
#> 2 chr2    200   400 +          2   170   200 0.85 
#> 3 chr2    300   500 -          2   130   200 0.65 
#> 4 chr2    800   900 -          0     0   100 0    
```
