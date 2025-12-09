# Create transcription start site features.

Create transcription start site features.

## Usage

``` r
create_tss(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md) in
  BED format

## See also

Other feature functions:
[`create_introns()`](https://rnabioco.github.io/valr/dev/reference/create_introns.md),
[`create_utrs3()`](https://rnabioco.github.io/valr/dev/reference/create_utrs3.md),
[`create_utrs5()`](https://rnabioco.github.io/valr/dev/reference/create_utrs5.md)

## Examples

``` r
x <- read_bed12(valr_example("mm9.refGene.bed.gz"))

create_tss(x)
#> # A tibble: 100 × 6
#>    chrom   start     end name         score strand
#>    <chr>   <dbl>   <dbl> <chr>        <chr> <chr> 
#>  1 chr1  3661578 3661579 NM_001011874 0     -     
#>  2 chr1  4399321 4399322 NM_001195662 0     -     
#>  3 chr1  4847774 4847775 NM_001159750 0     +     
#>  4 chr1  4847774 4847775 NM_011541    0     +     
#>  5 chr1  4848408 4848409 NM_001159751 0     +     
#>  6 chr1  5008814 5008815 NM_001290372 0     -     
#>  7 chr1  5009619 5009620 NM_021374    0     -     
#>  8 chr1  5060365 5060366 NM_001177795 0     -     
#>  9 chr1  5073166 5073167 NM_001310442 0     +     
#> 10 chr1  5073166 5073167 NM_133826    0     +     
#> # ℹ 90 more rows
```
