# Identify intersecting intervals.

Report intersecting intervals from `x` and `y` tbls.

## Usage

``` r
bed_intersect(x, ..., invert = FALSE, suffix = c(".x", ".y"), min_overlap = 1L)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- ...:

  one or more (e.g. a list of) `y`
  [`ivl_df()`](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)s.
  Each `y` may also be a path or URL to a bigWig (`.bw`) or bigBed
  (`.bb`) file, in which case only the regions spanned by `x` are read
  from it (local files and `http(s)://` URLs are both supported),
  avoiding the cost of loading the entire file.

- invert:

  report `x` intervals not in `y`

- suffix:

  colname suffixes in output

- min_overlap:

  minimum overlap in base pairs required for the operation. Defaults to
  `1`, which excludes book-ended intervals (those that touch but do not
  overlap), matching bedtools behavior. Set to `0` to include book-ended
  intervals (the legacy valr behavior).

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) with
original columns from `x` and `y` suffixed with `.x` and `.y`, and a new
`.overlap` column with the extent of overlap for the intersecting
intervals.

If multiple `y` tbls are supplied, the `.source` contains variable names
associated with each interval. All original columns from the `y` are
suffixed with `.y` in the output.

If `...` contains named inputs (i.e `a = y, b = z` or
`list(a = y, b = z)`), then `.source` will contain supplied names (see
examples).

## Details

input tbls are grouped by `chrom` by default, and additional groups can
be added using
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
For example, grouping by `strand` will constrain analyses to the same
strand. To compare opposing strands across two tbls, strands on the `y`
tbl can first be inverted using
[`flip_strands()`](https://rnabioco.github.io/valr/dev/reference/flip_strands.md).

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>

Other multiple set operations:
[`bed_closest()`](https://rnabioco.github.io/valr/dev/reference/bed_closest.md),
[`bed_coverage()`](https://rnabioco.github.io/valr/dev/reference/bed_coverage.md),
[`bed_map()`](https://rnabioco.github.io/valr/dev/reference/bed_map.md),
[`bed_subtract()`](https://rnabioco.github.io/valr/dev/reference/bed_subtract.md),
[`bed_window()`](https://rnabioco.github.io/valr/dev/reference/bed_window.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 25, 50,
  "chr1", 100, 125
)

y <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 30,     75
)

bed_glyph(bed_intersect(x, y))


bed_glyph(bed_intersect(x, y, invert = TRUE))


x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 100,    500,
  "chr2", 200,    400,
  "chr2", 300,    500,
  "chr2", 800,    900
)

y <- tibble::tribble(
  ~chrom, ~start, ~end, ~value,
  "chr1", 150,    400,  100,
  "chr1", 500,    550,  100,
  "chr2", 230,    430,  200,
  "chr2", 350,    430,  300
)

bed_intersect(x, y)
#> # A tibble: 5 × 7
#>   chrom start.x end.x start.y end.y value.y .overlap
#>   <chr>   <dbl> <dbl>   <dbl> <dbl>   <dbl>    <int>
#> 1 chr1      100   500     150   400     100      250
#> 2 chr2      200   400     230   430     200      170
#> 3 chr2      200   400     350   430     300       50
#> 4 chr2      300   500     230   430     200      130
#> 5 chr2      300   500     350   430     300       80

bed_intersect(x, y, invert = TRUE)
#> # A tibble: 1 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr2    800   900

# start and end of each overlapping interval
res <- bed_intersect(x, y)
dplyr::mutate(res,
  start = pmax(start.x, start.y),
  end = pmin(end.x, end.y)
)
#> # A tibble: 5 × 9
#>   chrom start.x end.x start.y end.y value.y .overlap start   end
#>   <chr>   <dbl> <dbl>   <dbl> <dbl>   <dbl>    <int> <dbl> <dbl>
#> 1 chr1      100   500     150   400     100      250   150   400
#> 2 chr2      200   400     230   430     200      170   230   400
#> 3 chr2      200   400     350   430     300       50   350   400
#> 4 chr2      300   500     230   430     200      130   300   430
#> 5 chr2      300   500     350   430     300       80   350   430

z <- tibble::tribble(
  ~chrom, ~start, ~end, ~value,
  "chr1", 150,    400,  100,
  "chr1", 500,    550,  100,
  "chr2", 230,    430,  200,
  "chr2", 750,    900,  400
)

bed_intersect(x, y, z)
#> # A tibble: 9 × 8
#>   chrom start.x end.x start.y end.y value.y .source .overlap
#>   <chr>   <dbl> <dbl>   <dbl> <dbl>   <dbl> <chr>      <int>
#> 1 chr1      100   500     150   400     100 y            250
#> 2 chr1      100   500     150   400     100 z            250
#> 3 chr2      200   400     230   430     200 y            170
#> 4 chr2      200   400     230   430     200 z            170
#> 5 chr2      200   400     350   430     300 y             50
#> 6 chr2      300   500     230   430     200 y            130
#> 7 chr2      300   500     230   430     200 z            130
#> 8 chr2      300   500     350   430     300 y             80
#> 9 chr2      800   900     750   900     400 z            100

bed_intersect(x, exons = y, introns = z)
#> # A tibble: 9 × 8
#>   chrom start.x end.x start.y end.y value.y .source .overlap
#>   <chr>   <dbl> <dbl>   <dbl> <dbl>   <dbl> <chr>      <int>
#> 1 chr1      100   500     150   400     100 exons        250
#> 2 chr1      100   500     150   400     100 introns      250
#> 3 chr2      200   400     230   430     200 exons        170
#> 4 chr2      200   400     230   430     200 introns      170
#> 5 chr2      200   400     350   430     300 exons         50
#> 6 chr2      300   500     230   430     200 exons        130
#> 7 chr2      300   500     230   430     200 introns      130
#> 8 chr2      300   500     350   430     300 exons         80
#> 9 chr2      800   900     750   900     400 introns      100

# a list of tbl_intervals can also be passed
bed_intersect(x, list(exons = y, introns = z))
#> # A tibble: 9 × 8
#>   chrom start.x end.x start.y end.y value.y .source .overlap
#>   <chr>   <dbl> <dbl>   <dbl> <dbl>   <dbl> <chr>      <int>
#> 1 chr1      100   500     150   400     100 exons        250
#> 2 chr1      100   500     150   400     100 introns      250
#> 3 chr2      200   400     230   430     200 exons        170
#> 4 chr2      200   400     230   430     200 introns      170
#> 5 chr2      200   400     350   430     300 exons         50
#> 6 chr2      300   500     230   430     200 exons        130
#> 7 chr2      300   500     230   430     200 introns      130
#> 8 chr2      300   500     350   430     300 exons         80
#> 9 chr2      800   900     750   900     400 introns      100

# a bigWig/bigBed file path or `http(s)://` URL can be used in place of a
# `y` tbl; only the regions spanned by `x` are read from the file
x <- tibble::tribble(
  ~chrom,  ~start,   ~end,
  "chr1",  4800000, 4830000,
  "chr10", 4850000, 4860000
)

bed_intersect(x, valr_example("test.bb"), min_overlap = 1L)
#> # A tibble: 2 × 15
#>   chrom start.x   end.x start.y   end.y name.y    score.y strand.y thickStart.y
#>   <chr>   <dbl>   <dbl>   <dbl>   <dbl> <chr>       <int> <chr>           <int>
#> 1 chr1  4800000 4830000 4797973 4836816 testgene        1 +             4797973
#> 2 chr10 4850000 4860000 4848118 4880877 diffchrom       1 +             4848118
#> # ℹ 6 more variables: thickEnd.y <int>, reserved.y <int>, blockCount.y <int>,
#> #   blockSizes.y <chr>, chromStarts.y <chr>, .overlap <int>
```
