# Select intervals bounded by a genome.

Used to remove out-of-bounds intervals, or trim interval coordinates
using a `genome`.

## Usage

``` r
bound_intervals(x, genome, trim = FALSE)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- genome:

  [genome_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

- trim:

  adjust coordinates for out-of-bounds intervals

## Value

[ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## See also

Other utilities:
[`bed12_to_exons()`](https://rnabioco.github.io/valr/reference/bed12_to_exons.md),
[`bed_makewindows()`](https://rnabioco.github.io/valr/reference/bed_makewindows.md),
[`flip_strands()`](https://rnabioco.github.io/valr/reference/flip_strands.md),
[`interval_spacing()`](https://rnabioco.github.io/valr/reference/interval_spacing.md)

## Examples

``` r
x <- tibble::tribble(
  ~chrom, ~start, ~end,
  "chr1", -100,   500,
  "chr1", 100,    1e9,
  "chr1", 500,    1000
)

genome <- read_genome(valr_example("hg19.chrom.sizes.gz"))

# out-of-bounds are removed by default ...
bound_intervals(x, genome)
#> # A tibble: 1 × 3
#>   chrom start   end
#>   <chr> <dbl> <dbl>
#> 1 chr1    500  1000

# ... or can be trimmed within the bounds of a genome
bound_intervals(x, genome, trim = TRUE)
#> # A tibble: 3 × 3
#>   chrom start       end
#>   <chr> <dbl>     <dbl>
#> 1 chr1      0       500
#> 2 chr1    100 249250621
#> 3 chr1    500      1000
```
