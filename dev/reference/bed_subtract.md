# Subtract two sets of intervals.

Subtract `y` intervals from `x` intervals.

## Usage

``` r
bed_subtract(x, y, any = FALSE, min_overlap = NULL)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- any:

  remove any `x` intervals that overlap `y`

- min_overlap:

  minimum overlap in base pairs required for subtraction. Set to `1` to
  exclude book-ended intervals (matching bedtools behavior), or `0` to
  include them (legacy valr behavior). The default will change from `0`
  to `1` in a future version.

## Details

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/dev/reference/flip_strands.md).

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html>

Other multiple set operations:
[`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md),
[`bed_coverage()`](https://rnabioco.github.io/valr/dev/reference/bed_coverage.md),
[`bed_intersect()`](https://rnabioco.github.io/valr/dev/reference/bed_intersect.md),
[`bed_map()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md),
[`bed_window()`](https://rnabioco.github.io/valr/dev/reference/bed_window.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1,      100
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 50,     75
)

bed_glyph(bed_subtract(x, y))
#> Warning: The `min_overlap` argument of `bed_subtract()` is deprecated as of valr 0.8.0.
#> ℹ The default will change from 0 (book-ended intervals overlap) to 1 (strict
#>   overlap) in a future version.
#> ℹ Set `min_overlap = 0L` to keep the legacy behavior, or `min_overlap = 1L` for
#>   bedtools-compatible behavior.


x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,    200,
  "chr1", 250,    400,
  "chr1", 500,    600,
  "chr1", 1000,   1200,
  "chr1", 1300,   1500
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 150,    175,
  "chr1", 510,    525,
  "chr1", 550,    575,
  "chr1", 900,    1050,
  "chr1", 1150,   1250,
  "chr1", 1299,   1501
)

bed_subtract(x, y)
#> # A tibble: 7 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    100   150
#> 2 chr1    175   200
#> 3 chr1    250   400
#> 4 chr1    500   510
#> 5 chr1    525   550
#> 6 chr1    575   600
#> 7 chr1   1050  1150

bed_subtract(x, y, any = TRUE)
#> # A tibble: 1 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    250   400
```
