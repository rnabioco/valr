# Fisher's test to measure overlap between two sets of intervals.

Calculate Fisher's test on number of intervals that are shared and
unique between two sets of `x` and `y` intervals.

## Usage

``` r
bed_fisher(x, y, genome)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## Value

[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## Details

Interval statistics can be used in combination with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and [`dplyr::do()`](https://dplyr.tidyverse.org/reference/do.html) to
calculate statistics for subsets of data. See
[`vignette('interval-stats')`](https://rnabioco.github.io/valr/articles/interval-stats.md)
for examples.

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/fisher.html>

Other interval statistics:
[`bed_absdist()`](https://rnabioco.github.io/valr/reference/bed_absdist.md),
[`bed_jaccard()`](https://rnabioco.github.io/valr/reference/bed_jaccard.md),
[`bed_projection()`](https://rnabioco.github.io/valr/reference/bed_projection.md),
[`bed_reldist()`](https://rnabioco.github.io/valr/reference/bed_reldist.md)

## Examples

``` r
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- bed_random(genome, n = 1e4, seed = 1010486)
y <- bed_random(genome, n = 1e4, seed = 9203911)

bed_fisher(x, y, genome)
#> # A tibble: 1 × 6
#>   estimate p.value conf.low conf.high method                         alternative
#>      <dbl>   <dbl>    <dbl>     <dbl> <chr>                          <chr>      
#> 1    0.945   0.707    0.722      1.22 Fisher's Exact Test for Count… two.sided  
```
