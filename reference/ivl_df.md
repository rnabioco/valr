# Bed-like data.frame requirements for valr functions

Required column names for interval dataframes are `chrom`, `start` and
`end`. Internally interval dataframes are validated using
`check_interval()`

Required column names for genome dataframes are `chrom` and `size`.
Internally genome dataframes are validated using `check_genome()`.

## Usage

``` r
check_interval(x)

check_genome(x)
```

## Arguments

- x:

  A `data.frame` or
  [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)

## Examples

``` r
# using tibble
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 1, 50,
  "chr1", 10, 75,
  "chr1", 100, 120
)

check_interval(x)
#> # A tibble: 3 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1      1    50
#> 2 chr1     10    75
#> 3 chr1    100   120

# using base R data.frame
x <- data.frame(
  chrom = "chr1",
  start = 0,
  end = 100,
  stringsAsFactors = FALSE
)

check_interval(x)
#> # A tibble: 1 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1      0   100

# example genome input

x <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 1e6
)

check_genome(x)
#> # A tibble: 1 × 2
#>   chrom    size
#>   <chr>   <dbl>
#> 1 chr1  1000000
```
