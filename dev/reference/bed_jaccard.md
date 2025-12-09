# Calculate the Jaccard statistic for two sets of intervals.

Quantifies the extent of overlap between to sets of intervals in terms
of base-pairs. Groups that are shared between input are used to
calculate the statistic for subsets of data.

## Usage

``` r
bed_jaccard(x, y)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

## Value

tibble with the following columns:

- `len_i` length of the intersection in base-pairs

- `len_u` length of the union in base-pairs

- `jaccard` value of jaccard statistic

- `n_int` number of intersecting intervals between `x` and `y`

If inputs are grouped, the return value will contain one set of values
per group.

## Details

The Jaccard statistic takes values of `[0,1]` and is measured as:

\$\$ J(x,y) = \frac{\mid x \bigcap y \mid} {\mid x \bigcup y \mid} =
\frac{\mid x \bigcap y \mid} {\mid x \mid + \mid y \mid - \mid x \bigcap
y \mid} \$\$

Interval statistics can be used in combination with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and [`dplyr::do()`](https://dplyr.tidyverse.org/reference/do.html) to
calculate statistics for subsets of data. See
[`vignette('interval-stats')`](https://rnabioco.github.io/valr/dev/articles/interval-stats.md)
for examples.

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html>

Other interval statistics:
[`bed_absdist()`](https://rnabioco.github.io/valr/dev/reference/bed_absdist.md),
[`bed_fisher()`](https://rnabioco.github.io/valr/dev/reference/bed_fisher.md),
[`bed_projection()`](https://rnabioco.github.io/valr/dev/reference/bed_projection.md),
[`bed_reldist()`](https://rnabioco.github.io/valr/dev/reference/bed_reldist.md)

## Examples

``` r
genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

x <- bed_random(genome, seed = 1010486)
y <- bed_random(genome, seed = 9203911)

bed_jaccard(x, y)
#> # A tibble: 1 × 4
#>       len_i      len_u jaccard      n
#>       <dbl>      <dbl>   <dbl>  <dbl>
#> 1 236195198 1708776065   0.160 399986

# calculate jaccard per chromosome
bed_jaccard(
  dplyr::group_by(x, chrom),
  dplyr::group_by(y, chrom)
)
#> # A tibble: 24 × 5
#>    chrom    len_i     len_u jaccard     n
#>    <chr>    <dbl>     <dbl>   <dbl> <dbl>
#>  1 chr1  18939736 137346996   0.160 32158
#>  2 chr10 10525496  75211235   0.163 17830
#>  3 chr11 10380469  74656102   0.161 17499
#>  4 chr12 10145996  73725426   0.160 17163
#>  5 chr13  8867978  63740684   0.162 14994
#>  6 chr14  8046352  59027220   0.158 13644
#>  7 chr15  7796590  56525195   0.160 13240
#>  8 chr16  6907153  49872823   0.161 11652
#>  9 chr17  6184210  44916521   0.160 10484
#> 10 chr18  6048908  43251434   0.163 10133
#> # ℹ 14 more rows
```
