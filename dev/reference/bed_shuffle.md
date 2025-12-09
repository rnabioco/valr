# Shuffle input intervals.

Shuffle input intervals.

## Usage

``` r
bed_shuffle(
  x,
  genome,
  incl = NULL,
  excl = NULL,
  max_tries = 1000,
  within = FALSE,
  seed = 0
)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- incl:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) of
  included intervals

- excl:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) of
  excluded intervals

- max_tries:

  maximum tries to identify a bounded interval

- within:

  shuffle within chromosomes

- seed:

  seed for reproducible intervals

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html>

Other randomizing operations:
[`bed_random()`](https://rnabioco.github.io/valr/dev/reference/bed_random.md)

## Examples

``` r
genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1e6,
  "chr2", 2e6,
  "chr3", 4e6
)

x <- bed_random(genome, seed = 1010486)

bed_shuffle(x, genome, seed = 9830491)
#> # A tibble: 1,000,000 × 3
#>    chrom   start     end
#>    <chr>   <int>   <int>
#>  1 chr2  1463822 1464822
#>  2 chr2   619967  620967
#>  3 chr2  1769865 1770865
#>  4 chr2   203953  204953
#>  5 chr3  2119387 2120387
#>  6 chr2  1216667 1217667
#>  7 chr3  2109652 2110652
#>  8 chr2   213473  214473
#>  9 chr1   154156  155156
#> 10 chr3  2201278 2202278
#> # ℹ 999,990 more rows
```
