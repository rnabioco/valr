# Calculate coverage across a genome

This function is useful for calculating interval coverage across an
entire genome.

## Usage

``` r
bed_genomecov(x, genome, zero_depth = FALSE)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- zero_depth:

  If TRUE, report intervals with zero depth. Zero depth intervals will
  be reported with respect to groups.

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) with
the an additional column:

- `.depth` depth of interval coverage

## Details

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/dev/reference/flip_strands.md).

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html>

Other single set operations:
[`bed_cluster()`](https://rnabioco.github.io/valr/dev/reference/bed_cluster.md),
[`bed_complement()`](https://rnabioco.github.io/valr/dev/reference/bed_complement.md),
[`bed_flank()`](https://rnabioco.github.io/valr/dev/reference/bed_flank.md),
[`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md),
[`bed_partition()`](https://rnabioco.github.io/valr/dev/reference/bed_partition.md),
[`bed_shift()`](https://rnabioco.github.io/valr/dev/reference/bed_shift.md),
[`bed_slop()`](https://rnabioco.github.io/valr/dev/reference/bed_slop.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 20, 70, "+",
  "chr1", 50, 100, "-",
  "chr1", 200, 250, "+",
  "chr1", 220, 250, "+"
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 500,
  "chr2", 1000
)

bed_genomecov(x, genome)
#> # A tibble: 5 × 4
#>   chrom start   end .depth
#>   <chr> <dbl> <dbl>  <int>
#> 1 chr1     20    50      1
#> 2 chr1     50    70      2
#> 3 chr1     70   100      1
#> 4 chr1    200   220      1
#> 5 chr1    220   250      2

bed_genomecov(dplyr::group_by(x, strand), genome)
#> # A tibble: 4 × 5
#>   chrom start   end strand .depth
#>   <chr> <dbl> <dbl> <chr>   <int>
#> 1 chr1     20    70 +           1
#> 2 chr1    200   220 +           1
#> 3 chr1    220   250 +           2
#> 4 chr1     50   100 -           1

bed_genomecov(dplyr::group_by(x, strand), genome, zero_depth = TRUE)
#> # A tibble: 11 × 5
#>    chrom start   end strand .depth
#>    <chr> <dbl> <dbl> <chr>   <int>
#>  1 chr1      0    20 +           0
#>  2 chr1      0    50 -           0
#>  3 chr1     20    70 +           1
#>  4 chr1     50   100 -           1
#>  5 chr1     70   200 +           0
#>  6 chr1    100   500 -           0
#>  7 chr1    200   220 +           1
#>  8 chr1    220   250 +           2
#>  9 chr1    250   500 +           0
#> 10 chr2      0  1000 +           0
#> 11 chr2      0  1000 -           0
```
