# Read genome files.

Genome files (UCSC "chromSize" files) contain chromosome name and size
information. These sizes are used by downstream functions to identify
computed intervals that have coordinates outside of the genome bounds.

## Usage

``` r
read_genome(path)
```

## Arguments

- path:

  containing chrom/contig names and sizes, one-pair-per-line,
  tab-delimited

## Value

[genome_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md),
sorted by `size`

## Note

URLs to genome files can also be used.

## See also

Other read functions:
[`read_bed()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md),
[`read_vcf()`](https://rnabioco.github.io/valr/dev/reference/read_vcf.md)

## Examples

``` r
read_genome(valr_example("hg19.chrom.sizes.gz"))
#> # A tibble: 25 × 2
#>    chrom      size
#>    <chr>     <dbl>
#>  1 chr1  249250621
#>  2 chr2  243199373
#>  3 chr3  198022430
#>  4 chr4  191154276
#>  5 chr5  180915260
#>  6 chr6  171115067
#>  7 chr7  159138663
#>  8 chrX  155270560
#>  9 chr8  146364022
#> 10 chr9  141213431
#> # ℹ 15 more rows

if (FALSE) { # \dontrun{
# `read_genome` accepts a URL
read_genome("https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes")
} # }
```
