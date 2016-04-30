[![Build Status](https://travis-ci.com/jayhesselberth/Rbedtools.svg?token=Q9WRSyqYnpS7KpFfTscp&branch=master)](https://travis-ci.com/jayhesselberth/Rbedtools) [![Coverage Status](https://img.shields.io/codecov/c/github/jayhesselberth/Rbedtools/master.svg)](https://codecov.io/github/jayhesselberth/Rbedtools?branch=master)

Rbedtools - Genome arithmetic in R
==================================

`Rbedtools` provides methods to do interval manipulations **within the R environment**, enabling fast explorative analysis of genome-scale data.

Installation
============

`Rbedtools` can be installed from github:

``` r
> devtools::install_github('jayhesselberth/Rbedtools')
```

Overview
========

`Rbedtools` enables the analysis of genome intervals in R. To keep the package lightweight, `Rbedtools` only supports BED+ and bedGraph formats. Key parts are implemented in `Rcpp` for speed.

The goal of `Rbedtools` is to enable easy analysis of genome-scale data sets **within R**. Moreover, `Rbedtools` makes use of new R libraries like `dplyr` and its pipe operator(`%>%`) for an expressive syntax that makes genome analysis fun. So a workflow like [this](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp3-plot-transcription-factor-occupancy-surrounding-the-transcription-start-site) becomes:

``` r
library(Rbedtools)
library(ggplot2)

tss_intervals <- read_bed('intervals.bed.gz')
chip_signal <- read_bedgraph('signal.bg.gz')

tss_intervals %>%
  bed_flank(size = 1000) %>%
  bed_makewindows(chip_signal, genome, win_size = 50) %>%
  bed_map(chip_signal, map_value = sum(value)) %>%
  group_by(win_id.y) %>%
  summarize(value_mean = mean(map_value), value_var = var(map_value)) %>%
  ggplot(aes(x = win_id.y, y = map_value)) + geom_point()
```

The main methods have the similar names to their `BEDtools` counterparts, so should be familiar to those users.

Vignette
========

See the vignettes for a full description of the package and examples (`browseVignettes(package = "Rbedtools")`)
