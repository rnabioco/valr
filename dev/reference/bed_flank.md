# Create flanking intervals from input intervals.

Create flanking intervals from input intervals.

## Usage

``` r
bed_flank(
  x,
  genome,
  both = 0,
  left = 0,
  right = 0,
  fraction = FALSE,
  strand = FALSE,
  trim = FALSE,
  ...
)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

- both:

  number of bases on both sizes

- left:

  number of bases on left side

- right:

  number of bases on right side

- fraction:

  define flanks based on fraction of interval length

- strand:

  define `left` and `right` based on strand

- trim:

  adjust coordinates for out-of-bounds intervals

- ...:

  extra arguments (not used)

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

## See also

<https://bedtools.readthedocs.io/en/latest/content/tools/flank.html>

Other single set operations:
[`bed_cluster()`](https://rnabioco.github.io/valr/dev/reference/bed_cluster.md),
[`bed_complement()`](https://rnabioco.github.io/valr/dev/reference/bed_complement.md),
[`bed_genomecov()`](https://rnabioco.github.io/valr/dev/reference/bed_genomecov.md),
[`bed_merge()`](https://rnabioco.github.io/valr/dev/reference/bed_merge.md),
[`bed_partition()`](https://rnabioco.github.io/valr/dev/reference/bed_partition.md),
[`bed_shift()`](https://rnabioco.github.io/valr/dev/reference/bed_shift.md),
[`bed_slop()`](https://rnabioco.github.io/valr/dev/reference/bed_slop.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", 25, 50,
  "chr1", 100, 125
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 130
)

bed_glyph(bed_flank(x, genome, both = 20))


x <- tibble::tribble(
  ~chrom, ~start, ~end, ~name, ~score, ~strand,
  "chr1", 500,    1000, ".",   ".",    "+",
  "chr1", 1000,   1500, ".",   ".",    "-"
)

genome <- tibble::tribble(
  ~chrom, ~size,
  "chr1", 5000
)

bed_flank(x, genome, left = 100)
#> # A tibble: 2 × 6
#>   chrom start   end name  score strand
#>   <chr> <dbl> <dbl> <chr> <chr> <chr> 
#> 1 chr1    400   500 .     .     +     
#> 2 chr1    900  1000 .     .     -     

bed_flank(x, genome, right = 100)
#> # A tibble: 2 × 6
#>   chrom start   end name  score strand
#>   <chr> <dbl> <dbl> <chr> <chr> <chr> 
#> 1 chr1   1000  1100 .     .     +     
#> 2 chr1   1500  1600 .     .     -     

bed_flank(x, genome, both = 100)
#> # A tibble: 4 × 6
#>   chrom start   end name  score strand
#>   <chr> <dbl> <dbl> <chr> <chr> <chr> 
#> 1 chr1    400   500 .     .     +     
#> 2 chr1   1000  1100 .     .     +     
#> 3 chr1    900  1000 .     .     -     
#> 4 chr1   1500  1600 .     .     -     

bed_flank(x, genome, both = 0.5, fraction = TRUE)
#> # A tibble: 4 × 6
#>   chrom start   end name  score strand
#>   <chr> <dbl> <dbl> <chr> <chr> <chr> 
#> 1 chr1    250   500 .     .     +     
#> 2 chr1   1000  1250 .     .     +     
#> 3 chr1    750  1000 .     .     -     
#> 4 chr1   1500  1750 .     .     -     
```
