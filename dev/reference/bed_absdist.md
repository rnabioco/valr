# Compute absolute distances between intervals.

Computes the absolute distance between the midpoint of each `x` interval
and the midpoints of each closest `y` interval.

## Usage

``` r
bed_absdist(x, y, genome)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) with
`.absdist` and `.absdist_scaled` columns.

## Details

Absolute distances are scaled by the inter-reference gap for the
chromosome as follows. For `Q` query points and `R` reference points on
a chromosome, scale the distance for each query point `i` to the closest
reference point by the inter-reference gap for each chromosome. If an
`x` interval has no matching `y` chromosome, `.absdist` is `NA`.

\$\$d_i(x,y) = min_k(\|q_i - r_k\|)\frac{R}{Length\\ of\\
chromosome}\$\$

Both absolute and scaled distances are reported as `.absdist` and
`.absdist_scaled`.

Interval statistics can be used in combination with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and
[`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
to calculate statistics for subsets of data. See
[`vignette('interval-stats')`](https://rnabioco.github.io/valr/dev/articles/interval-stats.md)
for examples.

## See also

<https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002529>

Other interval statistics:
[`bed_fisher()`](https://rnabioco.github.io/valr/dev/reference/bed_fisher.md),
[`bed_jaccard()`](https://rnabioco.github.io/valr/dev/reference/bed_jaccard.md),
[`bed_projection()`](https://rnabioco.github.io/valr/dev/reference/bed_projection.md),
[`bed_reldist()`](https://rnabioco.github.io/valr/dev/reference/bed_reldist.md)

## Examples

``` r
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- bed_random(genome, seed = 1010486)
y <- bed_random(genome, seed = 9203911)

bed_absdist(x, y, genome)
#> # A tibble: 1,000,000 × 5
#>    chrom start   end .absdist .absdist_scaled
#>    <chr> <dbl> <dbl>    <dbl>           <dbl>
#>  1 chr1   5184  6184     1392           0.448
#>  2 chr1   7663  8663     1087           0.350
#>  3 chr1   9858 10858     1526           0.491
#>  4 chr1  13805 14805     2421           0.779
#>  5 chr1  14081 15081     2697           0.868
#>  6 chr1  16398 17398     1700           0.547
#>  7 chr1  17486 18486      612           0.197
#>  8 chr1  22063 23063      466           0.150
#>  9 chr1  22494 23494      897           0.289
#> 10 chr1  29351 30351     1143           0.368
#> # ℹ 999,990 more rows
```
