# Create intron features.

Numbers in the `score` column are intron numbers from 5' to 3'
independent of strand. I.e., the first introns for `+` and `-` strand
genes both have `score` values of `1`.

## Usage

``` r
create_introns(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md) in BED12
  format

## See also

Other feature functions:
[`create_tss()`](https://rnabioco.github.io/valr/reference/create_tss.md),
[`create_utrs3()`](https://rnabioco.github.io/valr/reference/create_utrs3.md),
[`create_utrs5()`](https://rnabioco.github.io/valr/reference/create_utrs5.md)

## Examples

``` r
x <- read_bed12(valr_example("mm9.refGene.bed.gz"))

create_introns(x)
#> # A tibble: 1,583 × 6
#>    chrom   start     end name         score strand
#>    <chr>   <dbl>   <dbl> <chr>        <dbl> <chr> 
#>  1 chr1  3207049 3411782 NM_001011874     2 -     
#>  2 chr1  3411982 3660632 NM_001011874     1 -     
#>  3 chr1  4283093 4341990 NM_001195662     3 -     
#>  4 chr1  4342162 4342282 NM_001195662     2 -     
#>  5 chr1  4342918 4399250 NM_001195662     1 -     
#>  6 chr1  4848057 4857550 NM_001159750     1 +     
#>  7 chr1  4848057 4857550 NM_011541        1 +     
#>  8 chr1  4848584 4857550 NM_001159751     1 +     
#>  9 chr1  4857613 4868107 NM_001159750     2 +     
#> 10 chr1  4857613 4868107 NM_011541        2 +     
#> # ℹ 1,573 more rows
```
