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


x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 30,     75,
  "chr1", 50,     90,
  "chr1", 91,     120
)

bed_glyph(bed_merge(x))


bed_glyph(bed_cluster(x), label = ".id")

```
