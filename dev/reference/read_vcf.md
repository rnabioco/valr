# Read a VCF file.

Read a VCF file.

## Usage

``` r
read_vcf(vcf)
```

## Arguments

- vcf:

  vcf filename

## Value

`data_frame`

## Note

return value has `chrom`, `start` and `end` columns. Interval lengths
are the size of the 'REF' field.

## See also

Other read functions:
[`read_bed()`](https://rnabioco.github.io/valr/dev/reference/read_bed.md),
[`read_genome()`](https://rnabioco.github.io/valr/dev/reference/read_genome.md)

## Examples

``` r
vcf_file <- valr_example("test.vcf.gz")
read_vcf(vcf_file)
#> # A tibble: 11 × 18
#>    CHROM   POS ID      REF   ALT    QUAL FILTER INFO    FORMAT X1    X2    X3   
#>    <chr> <dbl> <chr>   <chr> <chr> <dbl> <chr>  <chr>   <chr>  <chr> <chr> <chr>
#>  1 1        10 1:10    A     T       100 PASS   NS=5;A… GT     0/0   ./.   0/0  
#>  2 1        20 1:20    G     C       100 PASS   NS=0;A… GT     ./.   ./.   ./.  
#>  3 1        30 1:30    C     A       100 PASS   NS=6;A… GT     0/0   0/0   0/0  
#>  4 1        40 1:40    A     C       100 PASS   NS=6;A… GT     0/0   0/0   0/0  
#>  5 1     10000 1:10000 G     C       100 PASS   NS=6;A… GT     0/0   0/0   0/0  
#>  6 1     20000 1:20000 T     A       100 PASS   NS=6;A… GT     1/1   1/1   1/1  
#>  7 4      5000 4:5000  A     T       100 PASS   NS=6;A… GT     1/1   1/1   1/1  
#>  8 4      6000 4:6000  C     T       100 PASS   NS=6;A… GT     1/1   1/1   1/1  
#>  9 X       800 X:800   A     C       100 PASS   NS=6;A… GT     1/1   1/1   1/1  
#> 10 X       900 X:900   A     T       100 PASS   NS=6;A… GT     1/1   1/1   1/1  
#> 11 X      1000 X:1000  T     G       100 PASS   NS=5;A… GT     1/1   1/1   1/1  
#> # ℹ 6 more variables: X4 <chr>, X5 <chr>, X6 <chr>, chrom <chr>, start <dbl>,
#> #   end <dbl>
```
