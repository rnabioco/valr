# Create 3' UTR features.

Create 3' UTR features.

## Usage

``` r
create_utrs3(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) in
  BED12 format

## See also

Other feature functions:
[`create_introns()`](https://rnabioco.github.io/valr/dev/reference/create_introns.md),
[`create_tss()`](https://rnabioco.github.io/valr/dev/reference/create_tss.md),
[`create_utrs5()`](https://rnabioco.github.io/valr/dev/reference/create_utrs5.md)

## Examples

``` r
x <- read_bed12(valr_example("mm9.refGene.bed.gz"))

create_utrs3(x)
#> # A tibble: 99 × 6
#>    chrom   start     end name         score strand
#>    <chr>   <dbl>   <dbl> <chr>        <chr> <chr> 
#>  1 chr1  3204562 3206102 NM_001011874 0     -     
#>  2 chr1  4280926 4283061 NM_001195662 0     -     
#>  3 chr1  4886445 4887990 NM_001159750 0     +     
#>  4 chr1  4886445 4887990 NM_011541    0     +     
#>  5 chr1  4886445 4887990 NM_001159751 0     +     
#>  6 chr1  4899656 4900554 NM_001290372 0     -     
#>  7 chr1  4899656 4900554 NM_021374    0     -     
#>  8 chr1  4899656 4900554 NM_001177795 0     -     
#>  9 chr1  5152246 5152630 NM_001310442 0     +     
#> 10 chr1  5152246 5152630 NM_133826    0     +     
#> # ℹ 89 more rows
```
