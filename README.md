read
====

``` r
library(Rbedtools)

bed_tbl <- read_bed('inst/extdata/3fields.bed.gz')
bed_tbl
#> Source: local data frame [10 x 3]
#> 
#>     chrom  start    end
#>    (fctr)  (int)  (int)
#> 1    chr1  11873  14409
#> 2    chr1  14361  19759
#> 3    chr1  14406  29370
#> 4    chr1  34610  36081
#> 5    chr1  69090  70008
#> 6    chr1 134772 140566
#> 7    chr1 321083 321115
#> 8    chr1 321145 321207
#> 9    chr1 322036 326938
#> 10   chr1 327545 328439

bedgraph_tbl <- read_bedgraph('inst/extdata/test.bg.gz')
bedgraph_tbl
#> Source: local data frame [4 x 4]
#> 
#>    chrom    start      end value
#>   (fctr)    (int)    (int) (dbl)
#> 1  chr19 49302000 49302300 -1.00
#> 2  chr19 49302300 49302600 -0.75
#> 3  chr19 49302600 49302900 -0.50
#> 4  chr19 49302900 49303200 -0.25
```

intersect
=========

``` r
# bed_tbl = "A" file
# bedgraph_tbl = "B" file
# A intersect B == B regions that intersect A intervals
bed_tbl %>%
  intersect(bedgraph_tbl) %>%
  summarize(n_insersects = n())
```
