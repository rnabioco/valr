# Calculate summaries from overlapping intervals.

Apply functions like [`min()`](https://rdrr.io/r/base/Extremes.html) and
[`max()`](https://rdrr.io/r/base/Extremes.html) to intersecting
intervals. `bed_map()` uses
[`bed_intersect()`](https://rnabioco.github.io/valr/reference/bed_intersect.md)
to identify intersecting intervals, so output columns will be suffixed
with `.x` and `.y`. Expressions that refer to input columns from `x` and
`y` columns must take these suffixes into account.

## Usage

``` r
bed_map(x, y, ..., min_overlap = 1L)

concat(.data, sep = ",")

values_unique(.data, sep = ",")

values(.data, sep = ",")
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- y:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- ...:

  name-value pairs specifying column names and expressions to apply

- min_overlap:

  minimum overlap in base pairs required for mapping. Default is `1`,
  meaning book-ended intervals (touching but not overlapping) are not
  included. Set to `0` to include book-ended intervals.

- .data:

  data

- sep:

  separator character

## Value

[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## Details

Non-intersecting intervals from `x` are included in the result with `NA`
values.

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/reference/flip_strands.md).

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/map.html>

Other multiple set operations:
[`bed_closest()`](https://rnabioco.github.io/valr/reference/bed_closest.md),
[`bed_coverage()`](https://rnabioco.github.io/valr/reference/bed_coverage.md),
[`bed_intersect()`](https://rnabioco.github.io/valr/reference/bed_intersect.md),
[`bed_subtract()`](https://rnabioco.github.io/valr/reference/bed_subtract.md),
[`bed_window()`](https://rnabioco.github.io/valr/reference/bed_window.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom ,
  ~start ,
  ~end   ,
  'chr1' ,
     100 ,
     250 ,
  'chr2' ,
     250 ,
     500
)

y <- tibble::tribble(
  ~chrom ,
  ~start ,
  ~end   ,
  ~value ,
  'chr1' ,
     100 ,
     250 ,
      10 ,
  'chr1' ,
     150 ,
     250 ,
      20 ,
  'chr2' ,
     250 ,
     500 ,
     500
)

bed_glyph(bed_map(x, y, value = sum(value)), label = 'value')


# summary examples
bed_map(x, y, .sum = sum(value))
#> # A tibble: 2 × 4
#>   chrom start   end  .sum
#>   <chr> <dbl> <dbl> <dbl>
#> 1 chr1    100   250    30
#> 2 chr2    250   500   500

bed_map(x, y, .min = min(value), .max = max(value))
#> # A tibble: 2 × 5
#>   chrom start   end  .min  .max
#>   <chr> <dbl> <dbl> <dbl> <dbl>
#> 1 chr1    100   250    10    20
#> 2 chr2    250   500   500   500

# identify non-intersecting intervals to include in the result
res <- bed_map(x, y, .sum = sum(value))
x_not <- bed_intersect(x, y, invert = TRUE)
dplyr::bind_rows(res, x_not)
#> # A tibble: 2 × 4
#>   chrom start   end  .sum
#>   <chr> <dbl> <dbl> <dbl>
#> 1 chr1    100   250    30
#> 2 chr2    250   500   500

# create a list-column
bed_map(x, y, .values = list(value))
#> # A tibble: 2 × 4
#>   chrom start   end .values  
#>   <chr> <dbl> <dbl> <list>   
#> 1 chr1    100   250 <dbl [2]>
#> 2 chr2    250   500 <dbl [1]>

# use `nth` family from dplyr
bed_map(x, y, .first = dplyr::first(value))
#> # A tibble: 2 × 4
#>   chrom start   end .first
#>   <chr> <dbl> <dbl>  <dbl>
#> 1 chr1    100   250     10
#> 2 chr2    250   500    500

bed_map(x, y, .absmax = abs(max(value)))
#> # A tibble: 2 × 4
#>   chrom start   end .absmax
#>   <chr> <dbl> <dbl>   <dbl>
#> 1 chr1    100   250      20
#> 2 chr2    250   500     500

bed_map(x, y, .count = length(value))
#> # A tibble: 2 × 4
#>   chrom start   end .count
#>   <chr> <dbl> <dbl>  <int>
#> 1 chr1    100   250      2
#> 2 chr2    250   500      1

bed_map(x, y, .vals = values(value))
#> # A tibble: 2 × 4
#>   chrom start   end .vals
#>   <chr> <dbl> <dbl> <chr>
#> 1 chr1    100   250 10,20
#> 2 chr2    250   500 500  

# count defaults are NA not 0; differs from bedtools2 ...
bed_map(x, y, .counts = dplyr::n())
#> # A tibble: 2 × 4
#>   chrom start   end .counts
#>   <chr> <dbl> <dbl>   <int>
#> 1 chr1    100   250       2
#> 2 chr2    250   500       1

# ... but NA counts can be coverted to 0's
dplyr::mutate(
  bed_map(x, y, .counts = dplyr::n()),
  .counts = ifelse(is.na(.counts), 0, .counts)
)
#> # A tibble: 2 × 4
#>   chrom start   end .counts
#>   <chr> <dbl> <dbl>   <int>
#> 1 chr1    100   250       2
#> 2 chr2    250   500       1
```
