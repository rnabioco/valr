# Convert BED12 to individual exons in BED6.

After conversion to BED6 format, the `score` column contains the exon
number, with respect to strand (i.e., the first exon for `-` strand
genes will have larger start and end coordinates).

## Usage

``` r
bed12_to_exons(x)
```

## Arguments

- x:

  [ivl_df](https://rnabioco.github.io/valr/reference/ivl_df.md)

## See also

Other utilities:
[`bed_makewindows()`](https://rnabioco.github.io/valr/reference/bed_makewindows.md),
[`bound_intervals()`](https://rnabioco.github.io/valr/reference/bound_intervals.md),
[`flip_strands()`](https://rnabioco.github.io/valr/reference/flip_strands.md),
[`interval_spacing()`](https://rnabioco.github.io/valr/reference/interval_spacing.md)

## Examples

``` r
x <- read_bed12(valr_example("mm9.refGene.bed.gz"))

bed12_to_exons(x)
#> # A tibble: 1,683 × 6
#>    chrom   start     end name         score strand
#>    <chr>   <dbl>   <dbl> <chr>        <int> <chr> 
#>  1 chr1  3204562 3207049 NM_001011874     3 -     
#>  2 chr1  3411782 3411982 NM_001011874     2 -     
#>  3 chr1  3660632 3661579 NM_001011874     1 -     
#>  4 chr1  4280926 4283093 NM_001195662     4 -     
#>  5 chr1  4341990 4342162 NM_001195662     3 -     
#>  6 chr1  4342282 4342918 NM_001195662     2 -     
#>  7 chr1  4399250 4399322 NM_001195662     1 -     
#>  8 chr1  4847774 4848057 NM_001159750     1 +     
#>  9 chr1  4847774 4848057 NM_011541        1 +     
#> 10 chr1  4848408 4848584 NM_001159751     1 +     
#> # ℹ 1,673 more rows
```
