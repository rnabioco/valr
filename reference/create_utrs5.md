# Create 5' UTR features.

Create 5' UTR features.

## Usage

``` r
create_utrs5(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md) in BED12
  format

## See also

Other feature functions:
[`create_introns()`](https://rnabioco.github.io/valr/reference/create_introns.md),
[`create_tss()`](https://rnabioco.github.io/valr/reference/create_tss.md),
[`create_utrs3()`](https://rnabioco.github.io/valr/reference/create_utrs3.md)

## Examples

``` r
x <- read_bed12(valr_example("mm9.refGene.bed.gz"))

create_utrs5(x)
#> # A tibble: 97 × 6
#>    chrom   start     end name         score strand
#>    <chr>   <dbl>   <dbl> <chr>        <chr> <chr> 
#>  1 chr1  3661429 3661579 NM_001011874 0     -     
#>  2 chr1  4399268 4399322 NM_001195662 0     -     
#>  3 chr1  4847774 4847994 NM_001159750 0     +     
#>  4 chr1  4847774 4847994 NM_011541    0     +     
#>  5 chr1  4848408 4848488 NM_001159751 0     +     
#>  6 chr1  4914046 5008815 NM_001290372 0     -     
#>  7 chr1  5009460 5009620 NM_021374    0     -     
#>  8 chr1  5060258 5060366 NM_001177795 0     -     
#>  9 chr1  5073166 5074531 NM_001310442 0     +     
#> 10 chr1  5073166 5074531 NM_133826    0     +     
#> # ℹ 87 more rows
```
