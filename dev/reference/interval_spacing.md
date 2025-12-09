# Calculate interval spacing.

Spacing for the first interval of each chromosome is undefined (`NA`).
The leading interval of an overlapping interval pair has a negative
value.

## Usage

``` r
interval_spacing(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) with
`.spacing` column.

## See also

Other utilities:
[`bed12_to_exons()`](https://rnabioco.github.io/valr/dev/reference/bed12_to_exons.md),
[`bed_makewindows()`](https://rnabioco.github.io/valr/dev/reference/bed_makewindows.md),
[`bound_intervals()`](https://rnabioco.github.io/valr/dev/reference/bound_intervals.md),
[`flip_strands()`](https://rnabioco.github.io/valr/dev/reference/flip_strands.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1,      100,
  "chr1", 150,    200,
  "chr2", 200,    300
)

interval_spacing(x)
#> # A tibble: 3 × 4
#>   chrom start   end .spacing
#>   <chr> <dbl> <dbl>    <dbl>
#> 1 chr1      1   100       NA
#> 2 chr1    150   200       50
#> 3 chr2    200   300       NA
```
