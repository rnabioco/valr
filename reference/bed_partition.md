# Partition intervals into elemental intervals

Convert a set of intervals into elemental intervals that contain each
start and end position in the set.

## Usage

``` r
bed_partition(x, ...)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- ...:

  name-value pairs specifying column names and expressions to apply

## Value

[`ivl_df()`](https://rnabioco.github.io/valr/reference/ivl_df.md)

## Details

Summary operations, such as
[`min()`](https://rdrr.io/r/base/Extremes.html) or
[`max()`](https://rdrr.io/r/base/Extremes.html) can be performed on
elemental intervals by specifying name-value pairs.

This function is useful for calculating summaries across overlapping
intervals without merging the intervals.

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/reference/flip_strands.md).

## See also

<https://bedops.readthedocs.io/en/latest/content/reference/set-operations/bedops.html#partition-p-partition>

Other single set operations:
[`bed_cluster()`](https://rnabioco.github.io/valr/reference/bed_cluster.md),
[`bed_complement()`](https://rnabioco.github.io/valr/reference/bed_complement.md),
[`bed_flank()`](https://rnabioco.github.io/valr/reference/bed_flank.md),
[`bed_genomecov()`](https://rnabioco.github.io/valr/reference/bed_genomecov.md),
[`bed_merge()`](https://rnabioco.github.io/valr/reference/bed_merge.md),
[`bed_shift()`](https://rnabioco.github.io/valr/reference/bed_shift.md),
[`bed_slop()`](https://rnabioco.github.io/valr/reference/bed_slop.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~value, ~strand,
  "chr1", 100, 500, 10, "+",
  "chr1", 200, 400, 20, "-",
  "chr1", 300, 550, 30, "+",
  "chr1", 550, 575, 2, "+",
  "chr1", 800, 900, 5, "+"
)


bed_glyph(bed_partition(x))

bed_glyph(bed_partition(x, value = sum(value)), label = "value")


bed_partition(x)
#> # A tibble: 7 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    100   200
#> 2 chr1    200   300
#> 3 chr1    300   400
#> 4 chr1    400   500
#> 5 chr1    500   550
#> 6 chr1    550   575
#> 7 chr1    800   900

# compute summary over each elemental interval
bed_partition(x, value = sum(value))
#> # A tibble: 7 × 4
#>   chrom start   end value
#>   <chr> <dbl> <dbl> <dbl>
#> 1 chr1    100   200    10
#> 2 chr1    200   300    30
#> 3 chr1    300   400    60
#> 4 chr1    400   500    40
#> 5 chr1    500   550    30
#> 6 chr1    550   575     2
#> 7 chr1    800   900     5

# partition and compute summaries based on group
x <- dplyr::group_by(x, strand)
bed_partition(x, value = sum(value))
#> # A tibble: 6 × 5
#>   chrom start   end strand value
#>   <chr> <dbl> <dbl> <chr>  <dbl>
#> 1 chr1    100   300 +         10
#> 2 chr1    200   400 -         20
#> 3 chr1    300   500 +         40
#> 4 chr1    500   550 +         30
#> 5 chr1    550   575 +          2
#> 6 chr1    800   900 +          5

# combine values across multiple tibbles
y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value, ~strand,
  "chr1", 10, 500, 100, "+",
  "chr1", 250, 420, 200, "-",
  "chr1", 350, 550, 300, "+",
  "chr1", 550, 555, 20, "+",
  "chr1", 800, 900, 50, "+"
)

x <- dplyr::bind_rows(x, y)
bed_partition(x, value = sum(value))
#> # A tibble: 11 × 5
#>    chrom start   end strand value
#>    <chr> <dbl> <dbl> <chr>  <dbl>
#>  1 chr1     10   100 +        100
#>  2 chr1    100   300 +        110
#>  3 chr1    200   250 -         20
#>  4 chr1    250   400 -        220
#>  5 chr1    300   350 +        140
#>  6 chr1    350   500 +        440
#>  7 chr1    400   420 -        200
#>  8 chr1    500   550 +        330
#>  9 chr1    550   555 +         22
#> 10 chr1    555   575 +          2
#> 11 chr1    800   900 +         55
```
