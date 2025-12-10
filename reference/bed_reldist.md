# Compute relative distances between intervals.

Compute relative distances between intervals.

## Usage

``` r
bed_reldist(x, y, detail = FALSE)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- detail:

  report relative distances for each `x` interval.

## Value

If `detail = FALSE`, a
[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md) that
summarizes calculated `.reldist` values with the following columns:

- `.reldist` relative distance metric

- `.counts` number of metric observations

- `.total` total observations

- `.freq` frequency of observation

If `detail = TRUE`, the `.reldist` column reports the relative distance
for each input `x` interval.

## Details

Interval statistics can be used in combination with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and [`dplyr::do()`](https://dplyr.tidyverse.org/reference/do.html) to
calculate statistics for subsets of data. See
[`vignette('interval-stats')`](https://rnabioco.github.io/valr/articles/interval-stats.md)
for examples.

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/reldist.html>

Other interval statistics:
[`bed_absdist()`](https://rnabioco.github.io/valr/reference/bed_absdist.md),
[`bed_fisher()`](https://rnabioco.github.io/valr/reference/bed_fisher.md),
[`bed_jaccard()`](https://rnabioco.github.io/valr/reference/bed_jaccard.md),
[`bed_projection()`](https://rnabioco.github.io/valr/reference/bed_projection.md)

## Examples

``` r
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- bed_random(genome, seed = 1010486)
y <- bed_random(genome, seed = 9203911)

bed_reldist(x, y)
#> # A tibble: 51 × 4
#>    .reldist .counts .total  .freq
#>       <dbl>   <int>  <int>  <dbl>
#>  1     0      20203 999939 0.0202
#>  2     0.01   20035 999939 0.0200
#>  3     0.02   19970 999939 0.0200
#>  4     0.03   19879 999939 0.0199
#>  5     0.04   20124 999939 0.0201
#>  6     0.05   20195 999939 0.0202
#>  7     0.06   20029 999939 0.0200
#>  8     0.07   20074 999939 0.0201
#>  9     0.08   20057 999939 0.0201
#> 10     0.09   20020 999939 0.0200
#> # ℹ 41 more rows

bed_reldist(x, y, detail = TRUE)
#> # A tibble: 999,939 × 4
#>    chrom start   end .reldist
#>    <chr> <dbl> <dbl>    <dbl>
#>  1 chr1   5184  6184   0.270 
#>  2 chr1   7663  8663   0.226 
#>  3 chr1   9858 10858   0.317 
#>  4 chr1  13805 14805   0.361 
#>  5 chr1  14081 15081   0.402 
#>  6 chr1  16398 17398   0.253 
#>  7 chr1  17486 18486   0.0912
#>  8 chr1  22063 23063   0.107 
#>  9 chr1  22494 23494   0.207 
#> 10 chr1  29351 30351   0.400 
#> # ℹ 999,929 more rows
```
