# Projection test for query interval overlap.

Projection test for query interval overlap.

## Usage

``` r
bed_projection(x, y, genome, by_chrom = FALSE)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- by_chrom:

  compute test per chromosome

## Value

[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md) with the
following columns:

- `chrom` the name of chromosome tested if `by_chrom = TRUE`, otherwise
  has a value of `whole_genome`

- `p.value` p-value from a binomial test. p-values \> 0.5 are converted
  to `1 - p-value` and `lower_tail` is `FALSE`

- `obs_exp_ratio` ratio of observed to expected overlap frequency

- `lower_tail` `TRUE` indicates the observed overlaps are in the lower
  tail of the distribution (e.g., less overlap than expected). `FALSE`
  indicates that the observed overlaps are in the upper tail of the
  distribution (e.g., more overlap than expected)

## Details

Interval statistics can be used in combination with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and [`dplyr::do()`](https://dplyr.tidyverse.org/reference/do.html) to
calculate statistics for subsets of data. See
[`vignette('interval-stats')`](https://rnabioco.github.io/valr/articles/interval-stats.md)
for examples.

## See also

<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529>

Other interval statistics:
[`bed_absdist()`](https://rnabioco.github.io/valr/reference/bed_absdist.md),
[`bed_fisher()`](https://rnabioco.github.io/valr/reference/bed_fisher.md),
[`bed_jaccard()`](https://rnabioco.github.io/valr/reference/bed_jaccard.md),
[`bed_reldist()`](https://rnabioco.github.io/valr/reference/bed_reldist.md)

## Examples

``` r
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- bed_random(genome, seed = 1010486)
y <- bed_random(genome, seed = 9203911)

bed_projection(x, y, genome)
#> # A tibble: 1 × 4
#>   chrom         p.value obs_exp_ratio lower_tail
#>   <chr>           <dbl>         <dbl> <chr>     
#> 1 whole_genome 0.000739          1.01 FALSE     

bed_projection(x, y, genome, by_chrom = TRUE)
#> # A tibble: 24 × 4
#>    chrom p.value obs_exp_ratio lower_tail
#>    <chr>   <dbl>         <dbl> <chr>     
#>  1 chr1   0.175          1.01  FALSE     
#>  2 chr10  0.0200         1.02  FALSE     
#>  3 chr11  0.111          1.01  FALSE     
#>  4 chr12  0.478          1.00  FALSE     
#>  5 chr13  0.241          1.01  FALSE     
#>  6 chr14  0.123          0.990 TRUE      
#>  7 chr15  0.338          1.00  FALSE     
#>  8 chr16  0.211          1.01  FALSE     
#>  9 chr17  0.171          0.990 TRUE      
#> 10 chr18  0.0272         1.02  FALSE     
#> # ℹ 14 more rows
```
