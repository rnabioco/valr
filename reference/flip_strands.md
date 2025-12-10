# Flip strands in intervals.

Flips positive (`+`) stranded intervals to negative (`-`) strands, and
vice-versa. Facilitates comparisons among intervals on opposing strands.

## Usage

``` r
flip_strands(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## See also

Other utilities:
[`bed12_to_exons()`](https://rnabioco.github.io/valr/reference/bed12_to_exons.md),
[`bed_makewindows()`](https://rnabioco.github.io/valr/reference/bed_makewindows.md),
[`bound_intervals()`](https://rnabioco.github.io/valr/reference/bound_intervals.md),
[`interval_spacing()`](https://rnabioco.github.io/valr/reference/interval_spacing.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end, ~strand,
  "chr1", 1,      100,  "+",
  "chr2", 1,      100,  "-"
)

flip_strands(x)
#> # A tibble: 2 × 4
#>   chrom start   end strand
#>   <chr> <dbl> <dbl> <chr> 
#> 1 chr1      1   100 -     
#> 2 chr2      1   100 +     
```
