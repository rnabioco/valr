# Cluster neighboring intervals.

The output `.id` column can be used in downstream grouping operations.
Default `max_dist = 0` means that both overlapping and book-ended
intervals will be clustered.

## Usage

``` r
bed_cluster(x, max_dist = 0)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- max_dist:

  maximum distance between clustered intervals.

## Value

[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md) with `.id`
column specifying sets of clustered intervals.

## Details

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/reference/flip_strands.md).

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/cluster.html>

Other single set operations:
[`bed_complement()`](https://rnabioco.github.io/valr/reference/bed_complement.md),
[`bed_flank()`](https://rnabioco.github.io/valr/reference/bed_flank.md),
[`bed_genomecov()`](https://rnabioco.github.io/valr/reference/bed_genomecov.md),
[`bed_merge()`](https://rnabioco.github.io/valr/reference/bed_merge.md),
[`bed_partition()`](https://rnabioco.github.io/valr/reference/bed_partition.md),
[`bed_shift()`](https://rnabioco.github.io/valr/reference/bed_shift.md),
[`bed_slop()`](https://rnabioco.github.io/valr/reference/bed_slop.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,    200,
  "chr1", 180,    250,
  "chr1", 250,    500,
  "chr1", 501,    1000,
  "chr2", 1,      100,
  "chr2", 150,    200
)

bed_cluster(x)
#> # A tibble: 6 × 4
#>   chrom start   end   .id
#>   <chr> <dbl> <dbl> <int>
#> 1 chr1    100   200     1
#> 2 chr1    180   250     1
#> 3 chr1    250   500     1
#> 4 chr1    501  1000     2
#> 5 chr2      1   100     3
#> 6 chr2    150   200     4

# glyph illustrating clustering of overlapping and book-ended intervals
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1,      10,
  "chr1", 5,      20,
  "chr1", 30,     40,
  "chr1", 40,     50,
  "chr1", 80,     90
)

bed_glyph(bed_cluster(x), label = ".id")

```
