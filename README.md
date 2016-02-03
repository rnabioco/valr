Installation
============

Rbedtools can be installed from github (CRAN submission is planned):

> devtools::install\_github('jayhesselberth/Rbedtools')

Usage and development
=====================

Problems? Submit an issue. Need a feature? Submit a pull request.

Overview
========

`Rbedtools` enables the analysis of genome intervals in R. To keep the package lightweight, `Rbedtools` only supports BED+ and bedGraph formats. Key parts are implemented in `Rcpp` for speed.

The goal of `Rbedtools` is to enable easy analysis of genome-scale data sets **within R**. So a workflow like [this](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp3-plot-transcription-factor-occupancy-surrounding-the-transcription-start-site) becomes:

``` r
tss_intervals %>%
  bed_flank(size = 1000) %>%
  bed_makewindows(size = 50) %>%
  bed_map(chip_signal) %>%
  ggplot(aes(x = win.num, y = signal), geom='point')
```

The main methods defined include:

-   `bed_intersect()`: intersections between sets of intervals
-   `bed_complement()`: identify intervals not covered by a query
-   `bed_merge()`: merge sets of intervals
-   `bed_map()`: analyze signals within intervals

-   `bed_sort()`: sort intervals. "If you make a BED file, sort the BED file."
-   `bed_makewindows()`: created labeled sub-intervals from a set of intervals
-   `bed12_to_exons()`: generate individual exons from BED12 format

Reading BED files
=================

Reading is done with `readr` for speed. Column types are coerced during reading.

Methods include:

-   `read_bed()`: read BED3 format
-   `read_bed12()`: read BED12 format
-   `read_bedgraph()`: read bedGraph format

``` r
library(Rbedtools)
#> Warning: replacing previous import by 'tidyr::%>%' when loading 'Rbedtools'

bed_path <- system.file('extdata/3fields.bed.gz', package = 'Rbedtools')
bed_tbl <- read_bed(bed_path)
bed_tbl
#> Source: local data frame [10 x 3]
#> 
#>     chrom  start    end
#>    (fctr)  (dbl)  (dbl)
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
```

Interval comparisons
====================

Interval intersections
----------------------

Interval intersections are implemented in `Rcpp` with the `chrom_sweep` algorithm from BEDtools.

``` r
# bed_tbl = "A" file
# bedgraph_tbl = "B" file
# A intersect B == B regions that intersect A intervals
bed_tbl %>%
  bed_intersect(bedgraph_tbl) %>%
  summarize(n_insersects = n())
```

Merging intervals
-----------------

Interval merging is implemented in `Rcpp` for speed.
