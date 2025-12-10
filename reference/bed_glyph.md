# Create example glyphs for valr functions.

Used to illustrate the output of valr functions with small examples.

## Usage

``` r
bed_glyph(expr, label = NULL)
```

## Arguments

- expr:

  expression to evaluate

- label:

  column name to use for label values. should be present in the result
  of the call.

## Value

[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 25,     50,
  "chr1", 100,    125
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value,
  "chr1", 30, 75, 50
)

bed_glyph(bed_intersect(x, y))
#> Warning: The `min_overlap` argument of `bed_intersect()` is deprecated as of valr 0.8.0.
#> ℹ The default will change from 0 (book-ended intervals overlap) to 1 (strict
#>   overlap) in a future version.
#> ℹ Set `min_overlap = 0L` to keep the legacy behavior, or `min_overlap = 1L` for
#>   bedtools-compatible behavior.


x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 30,     75,
  "chr1", 50,     90,
  "chr1", 91,     120
)

bed_glyph(bed_merge(x))


bed_glyph(bed_cluster(x), label = ".id")

```
