# Read BED and related files.

read functions for BED and related formats. Filenames can be local file
or URLs. The read functions load data into tbls with consistent `chrom`,
`start` and `end` colnames.

## Usage

``` r
read_bed(
  filename,
  col_types = bed12_coltypes,
  sort = TRUE,
  ...,
  n_fields = NULL
)

read_bed12(filename, ...)

read_bedgraph(filename, ...)

read_narrowpeak(filename, ...)

read_broadpeak(filename, ...)
```

## Arguments

- filename:

  file or URL

- col_types:

  column type spec for
  [`readr::read_tsv()`](https://readr.tidyverse.org/reference/read_delim.html)

- sort:

  sort the tbl by chrom and start

- ...:

  options to pass to
  [`readr::read_tsv()`](https://readr.tidyverse.org/reference/read_delim.html)

- n_fields:

  **\[deprecated\]**

## Value

[ivl_df](https://rnabioco.github.io/valr/dev/reference/ivl_df.md)

## Details

<https://genome.ucsc.edu/FAQ/FAQformat.html#format1>

<https://genome.ucsc.edu/FAQ/FAQformat.html#format1>

<https://genome.ucsc.edu/goldenPath/help/bedgraph.html>

<https://genome.ucsc.edu/FAQ/FAQformat.html#format12>

<https://genome.ucsc.edu/FAQ/FAQformat.html#format13>

## See also

Other read functions:
[`read_genome()`](https://rnabioco.github.io/valr/dev/reference/read_genome.md),
[`read_vcf()`](https://rnabioco.github.io/valr/dev/reference/read_vcf.md)

## Examples

``` r
# read_bed assumes 3 field BED format.
read_bed(valr_example("3fields.bed.gz"))
#> # A tibble: 10 × 3
#>    chrom  start    end
#>    <chr>  <dbl>  <dbl>
#>  1 chr1   11873  14409
#>  2 chr1   14361  19759
#>  3 chr1   14406  29370
#>  4 chr1   34610  36081
#>  5 chr1   69090  70008
#>  6 chr1  134772 140566
#>  7 chr1  321083 321115
#>  8 chr1  321145 321207
#>  9 chr1  322036 326938
#> 10 chr1  327545 328439

# result is sorted by chrom and start unless `sort = FALSE`
read_bed(valr_example("3fields.bed.gz"), sort = FALSE)
#> # A tibble: 10 × 3
#>    chrom  start    end
#>    <chr>  <int>  <int>
#>  1 chr1   11873  14409
#>  2 chr1   14361  19759
#>  3 chr1   14406  29370
#>  4 chr1   34610  36081
#>  5 chr1   69090  70008
#>  6 chr1  134772 140566
#>  7 chr1  321083 321115
#>  8 chr1  321145 321207
#>  9 chr1  322036 326938
#> 10 chr1  327545 328439


read_bed12(valr_example("mm9.refGene.bed.gz"))
#> # A tibble: 100 × 12
#>    chrom   start    end name  score strand cds_start cds_end item_rgb exon_count
#>    <chr>   <dbl>  <dbl> <chr> <chr> <chr>      <int>   <int> <chr>         <int>
#>  1 chr1  3204562 3.66e6 NM_0… 0     -        3206102 3661429 0                 3
#>  2 chr1  4280926 4.40e6 NM_0… 0     -        4283061 4399268 0                 4
#>  3 chr1  4847774 4.89e6 NM_0… 0     +        4847994 4886445 0                10
#>  4 chr1  4847774 4.89e6 NM_0… 0     +        4847994 4886445 0                10
#>  5 chr1  4848408 4.89e6 NM_0… 0     +        4848488 4886445 0                10
#>  6 chr1  4899656 5.01e6 NM_0… 0     -        4900554 4914046 0                 5
#>  7 chr1  4899656 5.01e6 NM_0… 0     -        4900554 5009460 0                 5
#>  8 chr1  4899656 5.06e6 NM_0… 0     -        4900554 5060258 0                 6
#>  9 chr1  5073166 5.15e6 NM_0… 0     +        5074531 5152246 0                13
#> 10 chr1  5073166 5.15e6 NM_1… 0     +        5074531 5152246 0                14
#> # ℹ 90 more rows
#> # ℹ 2 more variables: exon_sizes <chr>, exon_starts <chr>


read_bedgraph(valr_example("test.bg.gz"))
#> # A tibble: 4 × 4
#>   chrom    start      end value
#>   <chr>    <int>    <int> <dbl>
#> 1 chr19 49302000 49302300 -1   
#> 2 chr19 49302300 49302600 -0.75
#> 3 chr19 49302600 49302900 -0.5 
#> 4 chr19 49302900 49303200 -0.25


read_narrowpeak(valr_example("sample.narrowPeak.gz"))
#> # A tibble: 570 × 10
#>    chrom    start      end name  score strand signal pvalue qvalue  peak
#>    <chr>    <int>    <int> <chr> <int> <chr>   <dbl>  <dbl>  <dbl> <int>
#>  1 chr22 17372940 17373090 .         0 .           4   4.63     -1    -1
#>  2 chr22 17392200 17392350 .         0 .           5   4.67     -1    -1
#>  3 chr22 17398400 17398550 .         0 .          10  11.6      -1    -1
#>  4 chr22 17539180 17539330 .         0 .          21  30.9      -1    -1
#>  5 chr22 17652440 17652590 .         0 .           6   5.35     -1    -1
#>  6 chr22 17652780 17652930 .         0 .          12  12.5      -1    -1
#>  7 chr22 17980800 17980950 .         0 .          12  12.6      -1    -1
#>  8 chr22 18038260 18038410 .         0 .          29  36.0      -1    -1
#>  9 chr22 18225280 18225430 .         0 .          21  25.0      -1    -1
#> 10 chr22 18268020 18268170 .         0 .          14  13.0      -1    -1
#> # ℹ 560 more rows


read_broadpeak(valr_example("sample.broadPeak.gz"))
#> # A tibble: 1,181 × 9
#>    chrom    start      end name  score strand signal pvalue qvalue
#>    <chr>    <int>    <int> <chr> <int> <chr>   <dbl>  <dbl>  <dbl>
#>  1 chr22 16847903 16848440 .       503 .       10.5     2.5     -1
#>  2 chr22 16849452 16851326 .       483 .        9.81   15.7     -1
#>  3 chr22 16849955 16850086 .      1000 .       32.4     4.2     -1
#>  4 chr22 16850694 16850924 .       831 .       22.5     4.9     -1
#>  5 chr22 16852964 16853782 .       499 .       10.4     6.2     -1
#>  6 chr22 16855065 16855803 .       477 .        9.58    3.9     -1
#>  7 chr22 16855944 16856974 .       491 .       10.1     8.5     -1
#>  8 chr22 16857425 16857958 .       505 .       10.6     2.5     -1
#>  9 chr22 16858284 16858824 .       549 .       12.2     4.5     -1
#> 10 chr22 16859972 16862024 .       404 .        6.89    9.3     -1
#> # ℹ 1,171 more rows
```
