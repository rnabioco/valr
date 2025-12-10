# Generate randomly placed intervals on a genome.

Generate randomly placed intervals on a genome.

## Usage

``` r
bed_random(genome, length = 1000, n = 1e+06, seed = 0, sorted = TRUE)
```

## Arguments

- genome:

  [genome_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- length:

  length of intervals

- n:

  number of intervals to generate

- seed:

  seed RNG for reproducible intervals

- sorted:

  return sorted output

## Value

[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## Details

Sorting can be suppressed with `sorted = FALSE`.

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/random.html>

Other randomizing operations:
[`bed_shuffle()`](https://rnabioco.github.io/valr/reference/bed_shuffle.md)

## Examples

``` r
genome <- tibble::tribble(
  ~chrom,  ~size,
  "chr1",  10000000,
  "chr2",  50000000,
  "chr3",  60000000,
  "chrX",  5000000
)

bed_random(genome, seed = 10104)
#> # A tibble: 1,000,000 × 3
#>    chrom start   end
#>    <chr> <dbl> <dbl>
#>  1 chr1    265  1265
#>  2 chr1    315  1315
#>  3 chr1    513  1513
#>  4 chr1    625  1625
#>  5 chr1    635  1635
#>  6 chr1    653  1653
#>  7 chr1    731  1731
#>  8 chr1    859  1859
#>  9 chr1   1024  2024
#> 10 chr1   1038  2038
#> # ℹ 999,990 more rows

# sorting can be suppressed
bed_random(genome, sorted = FALSE, seed = 10104)
#> # A tibble: 1,000,000 × 3
#>    chrom    start      end
#>    <chr>    <dbl>    <dbl>
#>  1 chr3   2681468  2682468
#>  2 chr3  25364020 25365020
#>  3 chr2  43134407 43135407
#>  4 chr2  36344257 36345257
#>  5 chr3  39019690 39020690
#>  6 chr3  59869387 59870387
#>  7 chr3  56220373 56221373
#>  8 chr3  57965913 57966913
#>  9 chr2  25303342 25304342
#> 10 chr3  16018594 16019594
#> # ℹ 999,990 more rows

# 500 random intervals of length 500
bed_random(genome, length = 500, n = 500, seed = 10104)
#> # A tibble: 500 × 3
#>    chrom   start     end
#>    <chr>   <dbl>   <dbl>
#>  1 chr1   379360  379860
#>  2 chr1   394770  395270
#>  3 chr1  1215880 1216380
#>  4 chr1  1339287 1339787
#>  5 chr1  2046513 2047013
#>  6 chr1  2156755 2157255
#>  7 chr1  2189109 2189609
#>  8 chr1  2221665 2222165
#>  9 chr1  2223456 2223956
#> 10 chr1  2253135 2253635
#> # ℹ 490 more rows
```
