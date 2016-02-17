Rbedtools - Genome arithmetic in R
==================================

`Rbedtools` aims to fill a gap in genome analysis by providing interval manipulations within the R environment, enabling fast explorative analysis of genome-scale data.

Installation
============

Rbedtools can be installed from github:

    > devtools::install_github('jayhesselberth/Rbedtools')

Usage and development
=====================

Problems? Submit an [issue](https://github.com/jayhesselberth/Rbedtools/issues). Need a feature? Submit a [pull request](https://github.com/jayhesselberth/Rbedtools/pulls).

Overview
========

`Rbedtools` enables the analysis of genome intervals in R. To keep the package lightweight, `Rbedtools` only supports BED+ and bedGraph formats. Key parts are implemented in `Rcpp` for speed.

The goal of `Rbedtools` is to enable easy analysis of genome-scale data sets **within R**. So a workflow like [this](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp3-plot-transcription-factor-occupancy-surrounding-the-transcription-start-site) becomes:

``` r
tss_intervals %>%
  bed_flank(size = 1000) %>%
  bed_makewindows(win_size = 50) %>%
  bed_map(chip_signal) %>%
  ggplot(aes(x = win.num, y = signal)) + geom_point()
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

-   `read_bed()` \# default is BED3
-   `read_bed12()`
-   `read_bedgraph()`
-   `read_genome()` \# read chromSize files

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
