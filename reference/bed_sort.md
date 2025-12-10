# Sort a set of intervals.

Sort a set of intervals.

## Usage

``` r
bed_sort(x, by_size = FALSE, by_chrom = FALSE, reverse = FALSE)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- by_size:

  sort by interval size

- by_chrom:

  sort within chromosome

- reverse:

  reverse sort order

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/sort.html>

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr8", 500,    1000,
  "chr8", 1000,   5000,
  "chr8", 100,    200,
  "chr1", 100,    300,
  "chr1", 100,    200
)

# sort by chrom and start
bed_sort(x)
#> # A tibble: 5 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    100   200
#> 2 chr1    100   300
#> 3 chr8    100   200
#> 4 chr8    500  1000
#> 5 chr8   1000  5000

# reverse sort order
bed_sort(x, reverse = TRUE)
#> # A tibble: 5 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    100   300
#> 2 chr1    100   200
#> 3 chr8   1000  5000
#> 4 chr8    500  1000
#> 5 chr8    100   200

# sort by interval size
bed_sort(x, by_size = TRUE)
#> # A tibble: 5 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr8    100   200
#> 2 chr1    100   200
#> 3 chr1    100   300
#> 4 chr8    500  1000
#> 5 chr8   1000  5000

# sort by decreasing interval size
bed_sort(x, by_size = TRUE, reverse = TRUE)
#> # A tibble: 5 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr8   1000  5000
#> 2 chr8    500  1000
#> 3 chr1    100   300
#> 4 chr8    100   200
#> 5 chr1    100   200

# sort by interval size within chrom
bed_sort(x, by_size = TRUE, by_chrom = TRUE)
#> # A tibble: 5 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    100   200
#> 2 chr1    100   300
#> 3 chr8    100   200
#> 4 chr8    500  1000
#> 5 chr8   1000  5000
```
